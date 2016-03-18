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
import edu.washington.gs.maccoss.encyclopedia.algorithms.RetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackground;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
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
			
		} else if (arguments.containsKey("-pecan")) {
			Pecanpie.main(args);
			
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("Encyclopedia Help");
			Logger.logLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.logLine("Required Parameters: ");
			Logger.logLine("\t-i\tinput .DIA or .MZML file");
			Logger.logLine("\t-l\tlibrary .BLIB or .ELIB file");
			Logger.logLine("Other Parameters: ");
			Logger.logLine("\t-irtdb\tlibrary .IRTDB file");
			Logger.logLine("\t-o\toutput report file (default: [input file].pecan.txt)");
			
			TreeMap<String, String> defaults=new TreeMap<String, String>(PecanParameterParser.getDefaultParameters());
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
			
			String irtFileName=arguments.get("-irtdb");
			File irtFile=irtFileName==null?null:new File(irtFileName);

			File outputFile;
			if (arguments.containsKey("-o")) {
				outputFile=new File(arguments.get("-o"));
			} else {
				outputFile=new File(diaFile.getAbsolutePath()+".percolator.txt");
			}

			File featureFile=new File(diaFile.getAbsolutePath()+".features.txt");

			SearchParameters parameters=SearchParameterParser.parseParameters(arguments);
			LibraryScoringFactory factory=new EncyclopediaOneScoringFactory(parameters);
			Logger.logLine("Encyclopedia version "+factory.getVersion());

			Logger.logLine("Parameters:");
			Logger.logLine(" -i "+diaFile.getAbsolutePath());
			Logger.logLine(" -l "+libraryFile.getAbsolutePath());
			Logger.logLine(" -irtdb "+irtFileName);
			Logger.logLine(" -t "+arguments.get("-t"));
			Logger.logLine(" -o "+outputFile.getAbsolutePath());
			Logger.logLine(parameters.toString());

			try {
				LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile, Optional.ofNullable(irtFile));
				EncyclopediaJobData job=new EncyclopediaJobData(diaFile, library, featureFile, outputFile, factory);
				runSearch(new EmptyProgressIndicator(), job);
			} catch (Exception e) {
				System.err.println("Encountered Fatal Error!");
				e.printStackTrace();
			}
		}
	}

	public static void runSearch(ProgressIndicator progress, EncyclopediaJobData job) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		File diaFile=job.getDiaFile();
		File featureFile=job.getFeatureFile();
		File outputFile=job.getOutputFile();
		LibraryScoringFactory taskFactory=job.getTaskFactory();
		
		Logger.logLine("Converting files...");
		progress.update("Converting files...", Float.MIN_VALUE);
		
		SearchParameters parameters=taskFactory.getParameters();
		StripeFileInterface stripefile=MzmlToDIAConverter.getFile(diaFile, parameters);
		if (parameters.isDDA()) {
			EncyclopediaDDA.runSearch(progress, job.getLibrary(), stripefile, featureFile, outputFile, taskFactory);
		} else {
			runSearch(progress, job.getLibrary(), stripefile, featureFile, outputFile, taskFactory);
		}
	}
		
	public static void runSearch(ProgressIndicator progress, LibraryInterface library, StripeFileInterface stripefile, File featureFile, File outputFile, LibraryScoringFactory taskFactory) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		long startTime=System.currentTimeMillis();
		SearchParameters parameters=taskFactory.getParameters();
		
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

		PeptideScoringResultsConsumer writeResultsConsumer=taskFactory.getResultsConsumer(featureFile, new LinkedBlockingQueue<PeptideScoringResult>());
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
			Logger.logLine("Processing "+range+" ("+dutyCycle+")");
			
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
				tasks.add(entry.getReverse(parameters));
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

		ArrayList<ScoredObject<String>> passingPeptides=percolatePeptides(progress, featureFile, outputFile, taskFactory, saveResultsConsumer);
		
		Logger.logLine("Finished analysis! "+writeResultsConsumer.getNumberProcessed()+" total peaks processed, "+passingPeptides.size()+" peaks identified at "+(parameters.getPercolatorThreshold()*100f)+"% FDR ("+(Math.round((System.currentTimeMillis()-startTime)/1000f/6f)/10f)+" minutes)");
		Logger.logLine(""); 
	}

	public static ArrayList<ScoredObject<String>> percolatePeptides(ProgressIndicator progress, File featureFile, File outputFile, LibraryScoringFactory taskFactory, SaveResultsConsumer saveResultsConsumer) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		SearchParameters parameters=taskFactory.getParameters();
		
		ArrayList<ScoredObject<String>> passingPeptides;
		try {
			progress.update("Running Percolator ("+(parameters.getPercolatorThreshold()*100f)+"%)");
			File percolatorResultFile=new File(outputFile.getAbsolutePath()+".first_round.txt");
			
			passingPeptides=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorLocation(), featureFile, percolatorResultFile, parameters.getPercolatorThreshold());
			Logger.logLine("First pass: "+passingPeptides.size()+" peaks identified at "+(parameters.getPercolatorThreshold()*100f)+"% FDR");
			
			ArrayList<PeptideScoringResult> data=saveResultsConsumer.getSavedResults();
			RetentionTimeFilter filter=getRescoringModel(passingPeptides, data, outputFile);
			
			PeptideScoringResultsConsumer rescoredResultsConsumer=taskFactory.getResultsConsumer(featureFile, new LinkedBlockingQueue<PeptideScoringResult>());
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
			passingPeptides=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorLocation(), featureFile, outputFile, parameters.getPercolatorThreshold());
			
			progress.update(passingPeptides.size()+" peptides identified at "+(parameters.getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
			return passingPeptides;
			
		} catch (EncyclopediaException e) {
			Logger.errorLine("Fatal Error: "+e.getMessage());
			progress.update("Fatal Error: "+e.getMessage(), -1.0f);
			throw e;
		}
	}

	public static RetentionTimeFilter getRescoringModel(ArrayList<ScoredObject<String>> passingPeptides, ArrayList<PeptideScoringResult> data, File imageFileSeed) {
		HashSet<String> passingSeqs=new HashSet<String>();
		for (ScoredObject<String> pass : passingPeptides) {
			passingSeqs.add(pass.y);
		}
		
		HashSet<XYPoint> rtSet=new HashSet<XYPoint>();
		
		for (PeptideScoringResult result : data) {
			if (result.getGoodStripes().size()>0) {
				String peptideModSeq=result.getEntry().getPeptideModSeq();
				if (passingSeqs.contains(peptideModSeq+"+"+result.getEntry().getPrecursorCharge())) {
					LibraryEntry entry=result.getEntry();
					float entryTime=entry.getScanStartTime();

					Pair<ScoredObject<Stripe>, float[]> first=result.getGoodStripes().get(0);
					XYPoint point=new XYPoint(entryTime, first.x.y.getScanStartTime()/60.0f);
					rtSet.add(point);
				}
			}
		}
		ArrayList<XYPoint> rts=new ArrayList<XYPoint>(rtSet);
		Logger.logLine("Generating retention time mapping using "+rts.size()+" points...");
		RetentionTimeFilter filter=new RetentionTimeFilter(rts);
		
		filter.plot(rts, Optional.ofNullable(imageFileSeed));
		return filter;
	}
}
