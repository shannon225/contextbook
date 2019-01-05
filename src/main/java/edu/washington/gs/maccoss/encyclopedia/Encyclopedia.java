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
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackground;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
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
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.SaveResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.TeeResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.FileLogRecorder;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.VersioningDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class Encyclopedia {
	public static final String TARGET_LIBRARY_TAG="-l";
	public static final String OUTPUT_RESULT_TAG="-o";
	public static final String INPUT_DIA_TAG="-i";
	public static final String BACKGROUND_FASTA_TAG="-f";
	
	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		
		if (arguments.size()==0) {
			SearchGUIMain.runGUI(ProgramType.EncyclopeDIA);
			
		} else if (arguments.containsKey("-browser")) {
			DIABrowser.main(args);
		
		} else if (arguments.containsKey("-libexport")) {
			SearchToBLIB.main(args);
			
		} else if (arguments.containsKey("-walnut")||arguments.containsKey("-pecan")) {
			Walnut.main(args);
			
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("EncyclopeDIA Help");
			Logger.timelessLogLine("EncyclopeDIA is a library search engine for DIA data.");
			Logger.timelessLogLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .DIA or .MZML file");
			Logger.timelessLogLine("\t-f\tprotein .FASTA database");
			Logger.timelessLogLine("\t-l\tlibrary .ELIB file");
			Logger.timelessLogLine("Other Programs: ");
			Logger.timelessLogLine("\t-browser\trun ELIB Browser (use -browser -h for ELIB Browser help)");
			Logger.timelessLogLine("\t-libexport\trun Library Export (use -libexport -h for Library Export help)");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-o\toutput report file (default: [input file]"+EncyclopediaJobData.OUTPUT_FILE_SUFFIX+")");
			
			TreeMap<String, String> defaults=new TreeMap<String, String>(SearchParameterParser.getDefaultParameters());
			int maxWidth=0;
			for (String key : defaults.keySet()) {
				if (key.length()>maxWidth) maxWidth=key.length();
			}
			for (Entry<String, String> entry : defaults.entrySet()) {
				Logger.timelessLogLine("\t"+General.formatCellToWidth(entry.getKey(), maxWidth)+" (default: "+entry.getValue()+")");
			}
			System.exit(1);
			
		} else if (arguments.containsKey("-v")||arguments.containsKey("-version")||arguments.containsKey("--version")) {
			Logger.logLine("EncyclopeDIA version "+PecanOneScoringFactory.version);
			System.exit(1);
			
		} else {
			VersioningDetector.checkVersionCLI(ProgramType.EncyclopeDIA);
			
			if (!arguments.containsKey(INPUT_DIA_TAG)||!arguments.containsKey(TARGET_LIBRARY_TAG)||!arguments.containsKey(BACKGROUND_FASTA_TAG)) {
				Logger.errorLine("You are required to specify an input file ("+INPUT_DIA_TAG+"), a library file ("+TARGET_LIBRARY_TAG+"), and a fasta file ("+BACKGROUND_FASTA_TAG+")");
				System.exit(1);
			}

			File diaFile=new File(arguments.get(INPUT_DIA_TAG));
			File libraryFile=new File(arguments.get(TARGET_LIBRARY_TAG));
			File fastaFile=new File(arguments.get(Encyclopedia.BACKGROUND_FASTA_TAG));

			File outputFile;
			if (arguments.containsKey(OUTPUT_RESULT_TAG)) {
				outputFile=new File(arguments.get(OUTPUT_RESULT_TAG));
			} else {
				outputFile=new File(diaFile.getAbsolutePath()+EncyclopediaJobData.OUTPUT_FILE_SUFFIX);
			}

			try {
				FileLogRecorder logRecorder=new FileLogRecorder(new File(outputFile.getAbsolutePath()+EncyclopediaJobData.LOG_FILE_SUFFIX));
				Logger.addRecorder(logRecorder);
	
				SearchParameters parameters=SearchParameterParser.parseParameters(arguments);
				LibraryScoringFactory factory=new EncyclopediaOneScoringFactory(parameters);
				
				Logger.logLine("EncyclopeDIA version "+factory.getVersion());
	
				Logger.logLine("Parameters:");
				Logger.logLine(" "+INPUT_DIA_TAG+" "+diaFile.getAbsolutePath());
				Logger.logLine(" "+TARGET_LIBRARY_TAG+" "+libraryFile.getAbsolutePath());
				Logger.logLine(" "+OUTPUT_RESULT_TAG+" "+outputFile.getAbsolutePath());
				Logger.logLine(parameters.toString());

				LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
				EncyclopediaJobData job=new EncyclopediaJobData(diaFile, fastaFile, library, outputFile, factory);
				runSearch(new EmptyProgressIndicator(), job);
			} catch (Exception e) {
				Logger.errorLine("Encountered Fatal Error!");
				Logger.errorException(e);
			} finally {
				Logger.close();
			}
		}
	}

	public static void runSearch(ProgressIndicator progress, EncyclopediaJobData job) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		if (job.getPercolatorFiles().hasDataAvailable()) {
			try {
				ArrayList<PercolatorPeptide> passingPeptidesFromTSV=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideOutputFile(), job.getParameters(), false).x;
				
				File elibFile=job.getResultLibrary();
				if (!elibFile.exists()) {
					progress.update("Writing elib result library...");
					Logger.logLine("Writing elib result library...");
					ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
					jobs.add(job);
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
		
	static void runSearch(ProgressIndicator progress, EncyclopediaJobData job, StripeFileInterface stripefile) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		LibraryScoringFactory taskFactory=job.getTaskFactory();
		SearchParameters parameters=taskFactory.getParameters();
		if (parameters.getRtWindowInMin()>0.0f) {
			runIterativeSearch(progress, job, stripefile);
		} else {
			runTop5Search(progress, job, stripefile);
		}
	}

	static void runTop5Search(ProgressIndicator progress, EncyclopediaJobData job, StripeFileInterface stripefile) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		long startTime=System.currentTimeMillis();
		LibraryScoringFactory taskFactory=job.getTaskFactory();
		SearchParameters parameters=taskFactory.getParameters();

		Logger.logLine("Calculating features...");
		SaveResultsConsumer saveResultsConsumer=generateFeatureFile(progress, job, stripefile, Optional.empty());

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
		
		SearchToBLIB.convertElib(progress, job, elibFile, parameters);
		
		progress.update("Found "+passingPeptides.size()+" peptides identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
		Logger.logLine("Finished analysis! "+passingPeptides.size()+" peptides identified at "+(parameters.getPercolatorThreshold()*100f)+"% FDR ("+(Math.round((System.currentTimeMillis()-startTime)/1000f/6f)/10f)+" minutes)");
		Logger.logLine(""); 
	}
	
	static void runIterativeSearch(ProgressIndicator progress, EncyclopediaJobData job, StripeFileInterface stripefile) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		long startTime=System.currentTimeMillis();
		LibraryScoringFactory taskFactory=job.getTaskFactory();
		SearchParameters parameters=taskFactory.getParameters();
	
		Logger.logLine("Calculating features...");
		SaveResultsConsumer saveResultsConsumer=generateFeatureFile(progress, job, stripefile, Optional.empty());
	
		Logger.logLine("Running Percolator...");
		Pair<ArrayList<PercolatorPeptide>, RetentionTimeAlignmentInterface> percolatorResults=percolatePeptides(progress, job, stripefile, saveResultsConsumer);

		if (parameters.getScoringBreadthType().runRecalibration()) {
			Logger.logLine("Recalculating features...");
			saveResultsConsumer=generateFeatureFile(progress, job, stripefile, Optional.ofNullable(percolatorResults.y));
			percolatorResults=percolatePeptides(progress, job, stripefile, saveResultsConsumer);
			percolatorResults=repercolatePeptides(progress, job, stripefile, saveResultsConsumer, percolatorResults.y);
		}
		
		ArrayList<PercolatorPeptide> passingPeptides=percolatorResults.x;
		
		Logger.logLine("Writing elib result library...");
		File elibFile=job.getResultLibrary();
		ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
		jobs.add(job);
		
		SearchToBLIB.convertElib(progress, job, elibFile, parameters);
		
		progress.update("Found "+passingPeptides.size()+" peptides identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
		Logger.logLine("Finished analysis! "+passingPeptides.size()+" peptides identified at "+(parameters.getPercolatorThreshold()*100f)+"% FDR ("+(Math.round((System.currentTimeMillis()-startTime)/1000f/6f)/10f)+" minutes)");
		Logger.logLine(""); 
	}

	static SaveResultsConsumer generateFeatureFile(ProgressIndicator progress, EncyclopediaJobData job, StripeFileInterface stripefile, Optional<RetentionTimeAlignmentInterface> rtAlignment) throws IOException, SQLException, DataFormatException, InterruptedException {

		LibraryScoringFactory taskFactory=job.getTaskFactory();
		SearchParameters parameters=taskFactory.getParameters();
		LibraryInterface library=job.getLibrary();
		File featureFile=job.getPercolatorFiles().getInputTSV();
		
		int cores=parameters.getNumberOfThreadsUsed();

		Logger.logLine("Processing precursors scans...");
		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));

		// get targeted ranges
		ArrayList<Range> ranges=new ArrayList<Range>();
		for (Range range : stripefile.getRanges().keySet()) {
			if (!parameters.useTargetWindowCenter()||range.contains(parameters.getTargetWindowCenter())) {
				ranges.add(range);
			}
		}
		Collections.sort(ranges);

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
			String baseMessage="Working on "+range+" m/z";
			float baseIncrement=1.0f/numberOfTasks;
			float baseProgress=(1.0f+rangesFinished)/numberOfTasks;
			progress.update(baseMessage, baseProgress);
			
			float dutyCycle=stripefile.getRanges().get(range);
			Logger.logLine("Processing "+range+" m/z, ("+dutyCycle+" second duty cycle)");
			
			ArrayList<FragmentScan> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			Collections.sort(stripes);

			// prepare executor for background
			ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("STRIPE_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
			LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
			ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

			int count=0;
			ArrayList<LibraryEntry> entries=library.getEntries(range, true, parameters.getAAConstants());
			LibraryBackgroundInterface background=new LibraryBackground(entries);
			PSMScorer scorer=taskFactory.getLibraryScorer(background);
			
			for (LibraryEntry entry : entries) {
				count++;
				ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
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
				
				ArrayList<FragmentScan> localStripes;
				if (rtAlignment.isPresent()&&parameters.getRtWindowInMin()>0.0f) {
					float rtTargetInMin=rtAlignment.get().getYValue(entry.getRetentionTime()/60);
					localStripes=getScanSubsetFromStripes((rtTargetInMin-parameters.getRtWindowInMin())*60f, (rtTargetInMin+parameters.getRtWindowInMin())*60f, stripes);
				} else {
					localStripes=stripes;
				}
				executor.submit(taskFactory.getScoringTask(scorer, tasks, localStripes, range, dutyCycle, precursors, resultsQueue));
			}
			
			executor.shutdown();
			while (!executor.isTerminated()) {
				Logger.logLine(workQueue.size()+" peptides remaining for "+range+"...");
				float finishedFraction=(count-workQueue.size())/(float)count;
				progress.update(baseMessage, baseProgress+baseIncrement*(0.2f+finishedFraction*0.8f));
				Thread.sleep(500);
			}
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
			
			rangesFinished++;
		}
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

	public static Pair<ArrayList<PercolatorPeptide>, RetentionTimeAlignmentInterface> percolatePeptides(ProgressIndicator progress, EncyclopediaJobData job, StripeFileInterface stripefile, SaveResultsConsumer saveResultsConsumer) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		SearchParameters parameters=job.getParameters();
		
		try {
			progress.update("Running Percolator ("+(parameters.getPercolatorThreshold()*100f)+"%)");
			
			Pair<ArrayList<PercolatorPeptide>, Float> pair=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorVersionNumber(), job.getPercolatorFiles(), parameters.getEffectivePercolatorThreshold(), parameters.getAAConstants());
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
	public static Pair<ArrayList<PercolatorPeptide>, RetentionTimeAlignmentInterface> repercolatePeptides(ProgressIndicator progress, EncyclopediaJobData job, StripeFileInterface stripefile, SaveResultsConsumer saveResultsConsumer, RetentionTimeAlignmentInterface filter) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		SearchParameters parameters=job.getParameters();
		
		try {
			ArrayList<PeptideScoringResult> data=saveResultsConsumer.getSavedResults();
			PeptideScoringResultsConsumer rescoredResultsConsumer=job.getTaskFactory().getResultsConsumer(job.getPercolatorFiles().getInputTSV(), new LinkedBlockingQueue<PeptideScoringResult>(), stripefile);
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
			Pair<ArrayList<PercolatorPeptide>, Float> pair=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorVersionNumber(), job.getPercolatorFiles(), parameters.getEffectivePercolatorThreshold(), parameters.getAAConstants());
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

	public static RetentionTimeAlignmentInterface getRescoringModel(ArrayList<PercolatorPeptide> passingPeptides, ArrayList<PeptideScoringResult> data, EncyclopediaJobData job, boolean finalPass) {
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
		filter.plot(rts, Optional.ofNullable(new File(job.getPercolatorFiles().getPeptideOutputFile().getAbsolutePath()+passTag)));
		return filter;
	}
}
