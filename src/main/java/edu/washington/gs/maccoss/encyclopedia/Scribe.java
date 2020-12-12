package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.TreeMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RTRTPoint;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackground;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.scribe.ScribeJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.scribe.ScribeScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.WindowData;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.SaveResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.TeeResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.FileLogRecorder;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.VersioningDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ThreadableTask;
import gnu.trove.list.array.TFloatArrayList;

public class Scribe {
	public static final String TARGET_LIBRARY_TAG="-l";
	public static final String OUTPUT_RESULT_TAG="-o";
	public static final String INPUT_DIA_TAG="-i";
	public static final String BACKGROUND_FASTA_TAG="-f";
	
	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		
		if (arguments.size()==0) {
			SearchGUIMain.runGUI(ProgramType.EncyclopeDIA);

		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("Scribe Help");
			Logger.timelessLogLine("Scribe is a library search engine for DDA data.");
			Logger.timelessLogLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .DIA or .MZML file");
			Logger.timelessLogLine("\t-f\tprotein .FASTA database");
			Logger.timelessLogLine("\t-l\tlibrary .DLIB or .ELIB file");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-o\toutput report file (default: [input file]"+ScribeJobData.OUTPUT_FILE_SUFFIX+")");
			
			TreeMap<String, String> defaults=new TreeMap<String, String>(SearchParameterParser.getDefaultParameters());
			int maxWidth=0;
			for (String key : defaults.keySet()) {
				if (key.length()>maxWidth) maxWidth=key.length();
			}
			for (Entry<String, String> entry : defaults.entrySet()) {
				Logger.timelessLogLine("\t"+General.formatCellToWidth(entry.getKey(), maxWidth)+" (default: "+entry.getValue()+")");
			}

			Logger.timelessLogLine("\t"+Encyclopedia.QUIET_MODE_ARG+"\tsuppress log output to stdout/stderr");

			System.exit(1);
			
		} else if (arguments.containsKey("-v")||arguments.containsKey("-version")||arguments.containsKey("--version")) {
			Logger.logLine("Scribe version "+ProgramType.getGlobalVersion().toString());
			System.exit(1);
			
		} else {
			VersioningDetector.checkVersionCLI(ProgramType.EncyclopeDIA);
			
			if (!arguments.containsKey(INPUT_DIA_TAG)||!arguments.containsKey(TARGET_LIBRARY_TAG)||!arguments.containsKey(BACKGROUND_FASTA_TAG)) {
				Logger.errorLine("You are required to specify an input file ("+INPUT_DIA_TAG+"), a library file ("+TARGET_LIBRARY_TAG+"), and a fasta file ("+BACKGROUND_FASTA_TAG+")");
				System.exit(1);
			}

			File diaFile=new File(arguments.get(INPUT_DIA_TAG));
			File libraryFile=new File(arguments.get(TARGET_LIBRARY_TAG));
			File fastaFile=new File(arguments.get(Scribe.BACKGROUND_FASTA_TAG));

			File outputFile;
			if (arguments.containsKey(OUTPUT_RESULT_TAG)) {
				outputFile=new File(arguments.get(OUTPUT_RESULT_TAG));
			} else {
				outputFile=new File(diaFile.getAbsolutePath()+ScribeJobData.OUTPUT_FILE_SUFFIX);
			}

			try {
				if (arguments.containsKey(Encyclopedia.QUIET_MODE_ARG)) {
					Logger.PRINT_TO_SCREEN = false;
				}
				FileLogRecorder logRecorder=new FileLogRecorder(new File(outputFile.getAbsolutePath()+ScribeJobData.LOG_FILE_SUFFIX));
				Logger.addRecorder(logRecorder);
	
				SearchParameters parameters=SearchParameterParser.parseParameters(arguments);
				ScribeScoringFactory factory=new ScribeScoringFactory(parameters);
				
				Logger.logLine("Scribe version "+ProgramType.getGlobalVersion().toString());
	
				Logger.logLine("Parameters:");
				Logger.logLine(" "+INPUT_DIA_TAG+" "+diaFile.getAbsolutePath());
				Logger.logLine(" "+TARGET_LIBRARY_TAG+" "+libraryFile.getAbsolutePath());
				Logger.logLine(" "+OUTPUT_RESULT_TAG+" "+outputFile.getAbsolutePath());
				Logger.logLine(parameters.toString());

				LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
				ScribeJobData job=new ScribeJobData(diaFile, fastaFile, library, outputFile, factory);
				runSearch(new EmptyProgressIndicator(), job);
			} catch (Exception e) {
				Logger.errorLine("Encountered Fatal Error!");
				Logger.errorException(e);
			} finally {
				Logger.close();
			}
		}
	}

	public static void runSearch(ProgressIndicator progress, ScribeJobData job) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		if (getPercolatorData(job).hasDataAvailable()) {
			try {
				ArrayList<PercolatorPeptide> passingPeptidesFromTSV=PercolatorReader.getPassingPeptidesFromTSV(getPercolatorData(job).getPeptideOutputFile(), job.getParameters(), false).x;
				
				File elibFile=job.getResultLibrary();
				if (!elibFile.exists()) {
					progress.update("Writing elib result library...");
					Logger.logLine("Writing elib result library...");
					ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
					jobs.add(job);
					if (true) System.exit(1); // FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME
					SearchToBLIB.convert(progress, jobs, elibFile, false, false);
				}
				Logger.logLine("Previously found "+passingPeptidesFromTSV.size()+" peptides identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR");
				progress.update("Previously found "+passingPeptidesFromTSV.size()+" peptides identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
				//progress.update("Previously found "+passingPeptidesFromTSV.size()+" peptides ("+ParsimonyProteinGrouper.groupProteins(passingPeptidesFromTSV).size()+" proteins) identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);

				return;
			} catch (Exception e) {
				// problem! so just continue on and overwrite old result
				Logger.logLine("Found unexpected exception trying to read old results: ");
				Logger.logException(e);
				Logger.logLine("Just going to go ahead and reprocess this file!");
			}
		}
		
		Logger.logLine("Converting files...");
		progress.update("Converting files...", Float.MIN_VALUE);

		final StripeFileInterface stripefile = job.getDiaFileReader();
		runSearch(progress, job, stripefile);
		stripefile.close();
	}
		
	static void runSearch(ProgressIndicator progress, ScribeJobData job, StripeFileInterface stripefile) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		long startTime=System.currentTimeMillis();
		LibraryScoringFactory taskFactory=job.getTaskFactory();
		SearchParameters parameters=taskFactory.getParameters();

		Logger.logLine("Calculating features...");
		SaveResultsConsumer saveResultsConsumer=generateFeatureFile(progress, job, stripefile);

		Logger.logLine("Running Percolator...");
		Pair<ArrayList<PercolatorPeptide>, RetentionTimeAlignmentInterface> percolatorResults=percolatePeptides(progress, job, stripefile, saveResultsConsumer);
		if (parameters.getScoringBreadthType().runRecalibration()) {
			percolatorResults=repercolatePeptides(progress, job, stripefile, saveResultsConsumer, percolatorResults.y);
		}
		ArrayList<PercolatorPeptide> passingPeptides=percolatorResults.x;
		
		Logger.logLine("Writing elib result library...");
		File elibFile=job.getResultLibrary();
		ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
		jobs.add(job);
		
		//SearchToBLIB.convertElib(progress, job, elibFile, parameters);
		
		progress.update("Found "+passingPeptides.size()+" peptides identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
		Logger.logLine("Finished analysis! "+passingPeptides.size()+" peptides identified at "+(parameters.getPercolatorThreshold()*100f)+"% FDR ("+(Math.round((System.currentTimeMillis()-startTime)/1000f/6f)/10f)+" minutes)");
		Logger.logLine(""); 
	}
	
	public static final int NUMBER_OF_ISOTOPES_ABOVE_MONOISOTOPIC=1;
	private static final int NUMBER_OF_SPECTRA_IN_BATCH=1000;
	private static final float PERCENT_CONCURRENT_BATCHES_ABOVE_THREADCOUNT=0.5f; // keep 50% jobs waiting ready in the queue

	static SaveResultsConsumer generateFeatureFile(ProgressIndicator progress, ScribeJobData job, StripeFileInterface stripefile) throws IOException, SQLException, DataFormatException, InterruptedException {

		final LibraryScoringFactory taskFactory=job.getTaskFactory();
		final SearchParameters parameters=taskFactory.getParameters();
		final LibraryInterface library=job.getLibrary();
		File featureFile=getPercolatorData(job).getInputTSV();

		Logger.logLine("Processing precursors scans...");
		final PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));

		// get targeted ranges
		TFloatArrayList precursorList=new TFloatArrayList();
		for (Entry<Range, WindowData> entry : stripefile.getRanges().entrySet()) {
			float targetMZ=entry.getKey().getMiddle();
			for (int i = 0; i < entry.getValue().getNumberOfMSMS(); i++) {
				precursorList.add(targetMZ);
			}
		}
		precursorList.sort();
		
		ArrayList<Range> ranges=new ArrayList<>();
		int index=0;
		float lastStop=0.0f;
		while (true) {
			float start=Math.nextUp(lastStop);
			if (index+NUMBER_OF_SPECTRA_IN_BATCH>precursorList.size()) {
				float maxStop=Math.nextUp(precursorList.get(precursorList.size()-1));
				ranges.add(new Range(start, maxStop));
				break;
			} else {
				// random is deterministic, but it forces reads to be not be synchronizes
				index=index+RandomGenerator.randomIndex(NUMBER_OF_SPECTRA_IN_BATCH, index);
				lastStop=precursorList.get(index);
				if (lastStop<start) continue;

				ranges.add(new Range(start, lastStop));
			}
		}
		Logger.logLine("Found "+ranges.size()+" total ranges");

		// prepare executor for background
		int cores=parameters.getNumberOfThreadsUsed();
		int numberOfConcurrentJobs=Math.round(cores*(1f+PERCENT_CONCURRENT_BATCHES_ABOVE_THREADCOUNT));
		Logger.logLine("Preparing to maintain at most "+numberOfConcurrentJobs+" jobs for "+cores+" threads...");
		ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("Scribe thread factory").setDaemon(true).build();
		BlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
		ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

		PeptideScoringResultsConsumer writeResultsConsumer=taskFactory.getResultsConsumer(featureFile, new LinkedBlockingQueue<PeptideScoringResult>(), stripefile);
		SaveResultsConsumer saveResultsConsumer=new SaveResultsConsumer(new LinkedBlockingQueue<PeptideScoringResult>());
		
		BlockingQueue<PeptideScoringResult> resultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
		TeeResultsConsumer teeConsumer=new TeeResultsConsumer(resultsQueue, writeResultsConsumer, saveResultsConsumer);
		Thread consumer1Thread=new Thread(teeConsumer);
		Thread consumer2Thread=new Thread(writeResultsConsumer);
		Thread consumer3Thread=new Thread(saveResultsConsumer);
		consumer1Thread.start();
		consumer2Thread.start();
		consumer3Thread.start();

		// get stripes
		int rangesFinished=0;
		float numberOfTasks=2.0f+ranges.size();
		for (Range range : ranges) {
			float baseProgress=(1.0f+rangesFinished)/numberOfTasks;
			
			ThreadableTask<Nothing> task=new ThreadableTask<Nothing>() {
				@Override
				public String getTaskName() {
					return range.toString();
				}
				@Override
				protected Nothing process() {
					try {
						int count=0;
						double widerStart=range.getStart()-parameters.getFragmentTolerance().getTolerance(range.getStart());
						// assumes some +1Hs, so straight number of neutrons above target
						double targetStop=range.getStop()+MassConstants.neutronMass*NUMBER_OF_ISOTOPES_ABOVE_MONOISOTOPIC;
						double widerStop=targetStop+parameters.getFragmentTolerance().getTolerance(targetStop);
						ArrayList<LibraryEntry> entries=library.getEntries(new Range(widerStart, widerStop), true, parameters.getAAConstants());
						if (entries.size()==0) return Nothing.NOTHING;
						
						ArrayList<FragmentScan> stripes=stripefile.getStripes(range, -Float.MAX_VALUE, Float.MAX_VALUE, true);
						Collections.sort(stripes);
						
						String baseMessage="Working on "+range+" m/z ("+stripes.size()+" MS/MS, "+entries.size()+" total entries). Jobs in waiting: "+workQueue.size();
						progress.update(baseMessage, baseProgress);
						Logger.logLine(baseMessage);
						
						LibraryBackgroundInterface background=new LibraryBackground(entries);
						PSMScorer scorer=taskFactory.getLibraryScorer(background);
			
						ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
						for (LibraryEntry entry : entries) {
							count++;
							tasks.add(entry);
							tasks.add(entry.getDecoy(parameters));
							
							float extraDecoys=parameters.getNumberOfExtraDecoyLibrariesSearched();
							while (extraDecoys>0.0f) {
								if (extraDecoys<1.0f) {
									// check percentage
									float test=RandomGenerator.random(count);
									if (test>extraDecoys) {
										break;
									}
								}
								extraDecoys=extraDecoys-1.0f;
								LibraryEntry shuffle=entry.getShuffle(parameters, Float.hashCode(extraDecoys), false);
								tasks.add(shuffle);
								tasks.add(shuffle.getDecoy(parameters));
							}
						}
						//executor.submit(taskFactory.getDDAScoringTask(scorer, tasks, stripes, precursors, resultsQueue));
						return taskFactory.getDDAScoringTask(scorer, tasks, stripes, precursors, resultsQueue).call();
					} catch (IOException ioe) {
					} catch (SQLException sqle) {
					} catch (DataFormatException dfe) {
					}
			 return Nothing.NOTHING;
				}
			};
			System.err.println("Launched "+range);
			executor.submit(task);
			
			rangesFinished++;
		}
		executor.shutdown();
		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
		
		resultsQueue.put(PeptideScoringResult.POISON_RESULT);

		consumer1Thread.join();
		consumer2Thread.join();
		consumer3Thread.join();
		teeConsumer.close();
		progress.update("Organizing results", (1.0f+rangesFinished)/numberOfTasks);
		Logger.logLine(writeResultsConsumer.getNumberProcessed()+" total peptides processed.");
		return saveResultsConsumer;
	}

	
	static ArrayList<FragmentScan> getScanSubsetFromStripes(float minRT, float maxRT, ArrayList<FragmentScan> allScansInStripe) {
		ArrayList<FragmentScan> subset=new ArrayList<FragmentScan>();
		for (FragmentScan scan : allScansInStripe) {
			if (scan.getScanStartTime()>=minRT&&scan.getScanStartTime()<=maxRT) {
				subset.add(scan);
			}
		}
		return subset;
	}

	public static Pair<ArrayList<PercolatorPeptide>, RetentionTimeAlignmentInterface> percolatePeptides(ProgressIndicator progress, ScribeJobData job, StripeFileInterface stripefile, SaveResultsConsumer saveResultsConsumer) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		SearchParameters parameters=job.getParameters();
		
		try {
			progress.update("Running Percolator ("+(parameters.getPercolatorThreshold()*100f)+"%)");
			
			Pair<ArrayList<PercolatorPeptide>, Float> pair=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorVersionNumber(), getPercolatorData(job), parameters.getEffectivePercolatorThreshold(), parameters.getAAConstants());
			ArrayList<PercolatorPeptide> passingPeptides=pair.x;
			Logger.logLine("First pass: "+passingPeptides.size()+" peptides identified at "+(parameters.getPercolatorThreshold()*100f)+"% FDR");
			
			if (!parameters.getScoringBreadthType().runRecalibration()) {
				return new Pair<ArrayList<PercolatorPeptide>, RetentionTimeAlignmentInterface>(passingPeptides, null);
			}
			
			ArrayList<PeptideScoringResult> data=saveResultsConsumer.getSavedResults();
			RetentionTimeAlignmentInterface filter=getRescoringModel(passingPeptides, data, job, false);
			
			return new Pair<ArrayList<PercolatorPeptide>, RetentionTimeAlignmentInterface>(passingPeptides, filter);
		} catch (EncyclopediaException e) {
			Logger.errorLine("Fatal Error: "+e.getMessage());
			Logger.errorLine("Sorry, not feeling well today! Try again tomorrow!");
			progress.update("Fatal Error: "+e.getMessage(), -1.0f);
			throw e;
		}
	}

	private static PercolatorExecutionData getPercolatorData(ScribeJobData job) {
		return job.getPercolatorFiles().getDDAVersion();
	}
	
	public static Pair<ArrayList<PercolatorPeptide>, RetentionTimeAlignmentInterface> repercolatePeptides(ProgressIndicator progress, ScribeJobData job, StripeFileInterface stripefile, SaveResultsConsumer saveResultsConsumer, RetentionTimeAlignmentInterface filter) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		SearchParameters parameters=job.getParameters();
		
		try {
			ArrayList<PeptideScoringResult> data=saveResultsConsumer.getSavedResults();
			PeptideScoringResultsConsumer rescoredResultsConsumer=job.getTaskFactory().getResultsConsumer(getPercolatorData(job).getInputTSV(), new LinkedBlockingQueue<PeptideScoringResult>(), stripefile);
			Thread finalWriteConsumerThread=new Thread(rescoredResultsConsumer);
			finalWriteConsumerThread.start();
			BlockingQueue<PeptideScoringResult> resultList=rescoredResultsConsumer.getResultsQueue();
			for (PeptideScoringResult result : data) {
				PeptideScoringResult rescore=result.rescore(filter);
				resultList.add(rescore);
			}
			resultList.add(PeptideScoringResult.POISON_RESULT);
			finalWriteConsumerThread.join();
			rescoredResultsConsumer.close();
	
			progress.update("Re-running Percolator ("+(parameters.getPercolatorThreshold()*100f)+"%)");
			Pair<ArrayList<PercolatorPeptide>, Float> pair=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorVersionNumber(), getPercolatorData(job), parameters.getEffectivePercolatorThreshold(), parameters.getAAConstants());
			ArrayList<PercolatorPeptide> passingPeptides=pair.x;
			filter=getRescoringModel(passingPeptides, data, job, true);
			
			progress.update(passingPeptides.size()+" peptides identified at "+(parameters.getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
			return new Pair<ArrayList<PercolatorPeptide>, RetentionTimeAlignmentInterface>(passingPeptides, filter);
			
		} catch (EncyclopediaException e) {
			Logger.errorLine("Fatal Error: "+e.getMessage());
			Logger.errorLine("Sorry, not feeling well today! Try again tomorrow!");
			progress.update("Fatal Error: "+e.getMessage(), -1.0f);
			throw e;
		}
	}

	public static RetentionTimeAlignmentInterface getRescoringModel(ArrayList<PercolatorPeptide> passingPeptides, ArrayList<PeptideScoringResult> data, ScribeJobData job, boolean finalPass) {
		HashSet<String> passingSeqs=new HashSet<String>();
		for (PercolatorPeptide pass : passingPeptides) {
			passingSeqs.add(PercolatorPeptide.getPeptideData(pass.getPsmID()));
		}
		
		HashSet<XYPoint> rtSet=new HashSet<XYPoint>();
		
		for (PeptideScoringResult result : data) {
			if (result.getGoodStripes().size()>0) {
				String peptideModSeq=result.getEntry().getPeptideModSeq();
				if (passingSeqs.contains(peptideModSeq+"+"+result.getEntry().getPrecursorCharge())) {
					LibraryEntry entry=result.getEntry();
					float entryTime=entry.getScanStartTime();

					Pair<ScoredObject<FragmentScan>, float[]> first=result.getGoodStripes().get(0);
					XYPoint point=new RTRTPoint(entryTime/60.0f, first.x.y.getScanStartTime()/60.0f, entry.isDecoy(), entry.getPeptideModSeq());
					rtSet.add(point);
				}
			}
		}
		ArrayList<XYPoint> rts=new ArrayList<XYPoint>(rtSet);
		Logger.logLine("Generating retention time mapping using "+rts.size()+" points...");
		RetentionTimeAlignmentInterface filter=RetentionTimeFilter.getFilter(rts);
		
		final String passTag=finalPass?".final":".first";
		filter.plot(rts, Optional.ofNullable(new File(getPercolatorData(job).getPeptideOutputFile().getAbsolutePath()+passTag)));
		return filter;
	}
}
