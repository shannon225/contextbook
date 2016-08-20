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
import edu.washington.gs.maccoss.encyclopedia.algorithms.ParsimonyProteinGrouper;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackground;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ReverseLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.SaveResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.TeeResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class Encyclopedia {
	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		if (arguments.size()==0) {
			SearchGUIMain.runGUI(false);
			
		} else if (arguments.containsKey("-browser")) {
			DIABrowser.main(args);
		
		} else if (arguments.containsKey("-libexport")) {
			SearchToBLIB.main(args);
			
		} else if (arguments.containsKey("-pecan")) {
			Pecanpie.main(args);
			
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("Encyclopedia Help");
			Logger.logLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.logLine("Required Parameters: ");
			Logger.logLine("\t-i\tinput .DIA or .MZML file");
			Logger.logLine("\t-l\tlibrary .ELIB file");
			Logger.logLine("Other Parameters: ");
			Logger.logLine("\t-pecan\trun Pecanpie (use -pecan -h for Pecan help)");
			Logger.logLine("\t-browser\trun ELIB Browser (use -browser -h for ELIB Browser help)");
			Logger.logLine("\t-libexport\trun Library Export (use -libexport -h for Library Export help)");
			Logger.logLine("\t-o\toutput report file (default: [input file]"+EncyclopediaJobData.OUTPUT_FILE_SUFFIX+")");
			
			TreeMap<String, String> defaults=new TreeMap<String, String>(SearchParameterParser.getDefaultParameters());
			for (Entry<String, String> entry : defaults.entrySet()) {
				Logger.logLine("\t"+entry.getKey()+"\t(default: "+entry.getValue()+")");
			}
			System.exit(1);
			
		} else if (arguments.containsKey("-v")||arguments.containsKey("-version")||arguments.containsKey("--version")) {
			Logger.logLine("Encyclopedia version "+PecanOneScoringFactory.version);
			System.exit(1);
			
		} else {
			if (!arguments.containsKey("-i")||!arguments.containsKey("-l")) {
				Logger.errorLine("You are required to specify an input file (-i) and a library file (-l)");
				System.exit(1);
			}

			File diaFile=new File(arguments.get("-i"));
			File libraryFile=new File(arguments.get("-l"));

			File outputFile;
			if (arguments.containsKey("-o")) {
				outputFile=new File(arguments.get("-o"));
			} else {
				outputFile=new File(diaFile.getAbsolutePath()+EncyclopediaJobData.OUTPUT_FILE_SUFFIX);
			}

			SearchParameters parameters=SearchParameterParser.parseParameters(arguments);
			LibraryScoringFactory factory=new EncyclopediaOneScoringFactory(parameters);
			Logger.logLine("Encyclopedia version "+factory.getVersion());

			Logger.logLine("Parameters:");
			Logger.logLine(" -i "+diaFile.getAbsolutePath());
			Logger.logLine(" -l "+libraryFile.getAbsolutePath());
			Logger.logLine(" -t "+arguments.get("-t"));
			Logger.logLine(" -o "+outputFile.getAbsolutePath());
			Logger.logLine(parameters.toString());

			try {
				LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
				EncyclopediaJobData job=new EncyclopediaJobData(diaFile, library, outputFile, factory);
				runSearch(new EmptyProgressIndicator(), job);
			} catch (Exception e) {
				System.err.println("Encountered Fatal Error!");
				e.printStackTrace();
			}
		}
	}

	public static void runSearch(ProgressIndicator progress, EncyclopediaJobData job) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		File outputFile=job.getOutputFile();
		if (outputFile.exists()&&outputFile.canRead()) {
			try {
				ArrayList<PercolatorPeptide> passingPeptidesFromTSV=PercolatorReader.getPassingPeptidesFromTSV(outputFile, job.getParameters().getEffectivePercolatorThreshold());
				ArrayList<ProteinGroup> proteins=ParsimonyProteinGrouper.groupProtein(passingPeptidesFromTSV);
				
				File elibFile=job.getResultLibrary();
				if (!elibFile.exists()) {
					progress.update("Writing elib result library...");
					Logger.logLine("Writing elib result library...");
					ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
					jobs.add(job);
					SearchToBLIB.convert(progress, jobs, elibFile, false, true);
				}
				progress.update("Previously found "+passingPeptidesFromTSV.size()+" peptides ("+proteins.size()+" proteins) identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);

				return;
			} catch (Exception e) {
				// problem! so just continue on and overwrite old result
			}
		}
		
		File diaFile=job.getDiaFile();
		
		Logger.logLine("Converting files...");
		progress.update("Converting files...", Float.MIN_VALUE);
		
		SearchParameters parameters=job.getParameters();
		StripeFileInterface stripefile=MzmlToDIAConverter.getFile(diaFile, parameters);
		if (parameters.isDDA()) {
			EncyclopediaDDA.runSearch(progress, job, stripefile);
		} else {
			runSearch(progress, job, stripefile);
		}
		stripefile.close();
	}
		
	public static void runSearch(ProgressIndicator progress, EncyclopediaJobData job, StripeFileInterface stripefile) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		long startTime=System.currentTimeMillis();
		LibraryScoringFactory taskFactory=job.getTaskFactory();
		SearchParameters parameters=taskFactory.getParameters();
		LibraryInterface library=job.getLibrary();
		File featureFile=job.getFeatureFile();
		
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

		PeptideScoringResultsConsumer writeResultsConsumer=taskFactory.getResultsConsumer(featureFile, new LinkedBlockingQueue<PeptideScoringResult>(), stripefile.getFile());
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
		boolean everyOtherDecoy=true;
		for (Range range : ranges) {
			String baseMessage="Working on "+range+" m/z";
			float baseIncrement=1.0f/numberOfTasks;
			float baseProgress=(1.0f+rangesFinished)/numberOfTasks;
			progress.update(baseMessage, baseProgress);
			
			float dutyCycle=stripefile.getRanges().get(range);
			Logger.logLine("Processing "+range+" m/z, ("+dutyCycle+" second duty cycle)");
			
			ArrayList<Stripe> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			Collections.sort(stripes);

			// prepare executor for background
			ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("SWATH_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
			LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
			ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

			int count=0;
			ArrayList<LibraryEntry> entries=library.getEntries(range, true);
			LibraryBackground background=new LibraryBackground(entries);
			PSMScorer scorer=taskFactory.getLibraryScorer(background);
			
			for (LibraryEntry entry : entries) {
				count++;
				ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
				tasks.add(entry);
				tasks.add(entry.getDecoy(parameters, false));
				if (parameters.isAddExtraDecoys()) {
					if (everyOtherDecoy) {
						ReverseLibraryEntry shuffle=entry.getDecoy(parameters, true);
						tasks.add(shuffle);
						tasks.add(shuffle.getDecoy(parameters, false));
					}
					everyOtherDecoy=!everyOtherDecoy;
				}
				executor.submit(taskFactory.getScoringTask(scorer, tasks, stripes, dutyCycle, precursors, resultsQueue));
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

		Logger.logLine("Running Percolator...");
		ArrayList<PercolatorPeptide> passingPeptides=percolatePeptides(progress, job, stripefile, saveResultsConsumer);
		
		Logger.logLine("Writing elib result library...");
		File elibFile=job.getResultLibrary();
		ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
		jobs.add(job);
		SearchToBLIB.convert(progress, jobs, elibFile, false, true);

		Logger.logLine("Grouping proteins...");
		ArrayList<ProteinGroup> proteins=ParsimonyProteinGrouper.groupProtein(passingPeptides);
		
		progress.update("Found "+passingPeptides.size()+" peptides ("+proteins.size()+" proteins) identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
		Logger.logLine("Finished analysis! "+writeResultsConsumer.getNumberProcessed()+" total peptides processed, "+passingPeptides.size()+" peptides ("+proteins.size()+" proteins) identified at "+(parameters.getPercolatorThreshold()*100f)+"% FDR ("+(Math.round((System.currentTimeMillis()-startTime)/1000f/6f)/10f)+" minutes)");
		Logger.logLine(""); 
	}

	public static ArrayList<PercolatorPeptide> percolatePeptides(ProgressIndicator progress, EncyclopediaJobData job, StripeFileInterface stripefile, SaveResultsConsumer saveResultsConsumer) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		SearchParameters parameters=job.getParameters();
		
		ArrayList<PercolatorPeptide> passingPeptides;
		try {
			progress.update("Running Percolator ("+(parameters.getPercolatorThreshold()*100f)+"%)");
			File percolatorResultFile=job.getFirstPassPercolator();
			
			passingPeptides=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorLocation(), job.getFeatureFile(), percolatorResultFile, parameters.getEffectivePercolatorThreshold());
			Logger.logLine("First pass: "+passingPeptides.size()+" peptides identified at "+(parameters.getPercolatorThreshold()*100f)+"% FDR");
			
			ArrayList<PeptideScoringResult> data=saveResultsConsumer.getSavedResults();
			RetentionTimeFilter filter=getRescoringModel(passingPeptides, data, job);
			
			PeptideScoringResultsConsumer rescoredResultsConsumer=job.getTaskFactory().getResultsConsumer(job.getFeatureFile(), new LinkedBlockingQueue<PeptideScoringResult>(), stripefile.getFile());
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
			passingPeptides=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorLocation(), job.getFeatureFile(), job.getOutputFile(), parameters.getEffectivePercolatorThreshold());
			
			progress.update(passingPeptides.size()+" peptides identified at "+(parameters.getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
			return passingPeptides;
			
		} catch (EncyclopediaException e) {
			Logger.errorLine("Fatal Error: "+e.getMessage());
			Logger.errorLine("Sorry, not feeling well today! Try again tomorrow!");
			progress.update("Fatal Error: "+e.getMessage(), -1.0f);
			throw e;
		}
	}

	public static RetentionTimeFilter getRescoringModel(ArrayList<PercolatorPeptide> passingPeptides, ArrayList<PeptideScoringResult> data, EncyclopediaJobData job) {
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

					Pair<ScoredObject<Stripe>, float[]> first=result.getGoodStripes().get(0);
					XYPoint point=new XYPoint(entryTime/60.0f, first.x.y.getScanStartTime()/60.0f);
					rtSet.add(point);
				}
			}
		}
		ArrayList<XYPoint> rts=new ArrayList<XYPoint>(rtSet);
		Logger.logLine("Generating retention time mapping using "+rts.size()+" points...");
		RetentionTimeFilter filter=new RetentionTimeFilter(rts);
		
		filter.plot(rts, Optional.ofNullable(job.getOutputFile()));
		return filter;
	}
}
