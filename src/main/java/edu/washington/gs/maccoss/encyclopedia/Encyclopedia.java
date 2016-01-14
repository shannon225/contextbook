package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
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
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import gnu.trove.set.hash.TDoubleHashSet;

public class Encyclopedia {
	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		if (arguments.size()==0) {
			PecanMain.main(args);
			
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("Encyclopedia Help");
			Logger.logLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.logLine("Required Parameters: ");
			Logger.logLine("\t-i\tinput .DIA or .MZML file");
			Logger.logLine("\t-l\tlibrary .ELIB file");
			Logger.logLine("Other Parameters: ");
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

			File outputFile;
			if (arguments.containsKey("-o")) {
				outputFile=new File(arguments.get("-o"));
			} else {
				outputFile=new File(diaFile.getAbsolutePath()+".pecan.txt");
			}

			File featureFile=new File(outputFile.getAbsolutePath()+".features.txt");

			SearchParameters parameters=SearchParameterParser.parseParameters(arguments);
			LibraryScoringFactory factory=new EncyclopediaOneScoringFactory(parameters, featureFile);
			Logger.logLine("Encyclopedia version "+factory.getVersion());

			Logger.logLine("Parameters:");
			Logger.logLine(" -i "+diaFile.getAbsolutePath());
			Logger.logLine(" -l "+libraryFile.getAbsolutePath());
			Logger.logLine(" -t "+arguments.get("-t"));
			Logger.logLine(" -o "+outputFile.getAbsolutePath());
			Logger.logLine(parameters.toString());

			try {
				runSearch(new EmptyProgressIndicator(), libraryFile, diaFile, featureFile, outputFile, factory);
			} catch (Exception e) {
				System.err.println("Encountered Fatal Error!");
				e.printStackTrace();
			}
		}
	}
	
	static void runSearch(ProgressIndicator progress, File libraryFile, File diaFile, File featureFile, File outputFile, LibraryScoringFactory taskFactory) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		Logger.logLine("Converting files...");
		progress.update("Converting files...", Float.MIN_VALUE);
		
		LibraryFile library=new LibraryFile();
		library.openFile(libraryFile);

		SearchParameters parameters=taskFactory.getParameters();
		StripeFileInterface stripefile=MzmlToDIAConverter.getFile(diaFile, parameters);
		runSearch(progress, library, stripefile, featureFile, outputFile, taskFactory);
	}
		
	static void runSearch(ProgressIndicator progress, LibraryFile library, StripeFileInterface stripefile, File featureFile, File outputFile, LibraryScoringFactory taskFactory) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		long startTime=System.currentTimeMillis();
		PSMScorer pecanScorer=taskFactory.getLibraryScorer();
		SearchParameters parameters=taskFactory.getParameters();
		
		int cores=Runtime.getRuntime().availableProcessors();


		Logger.logLine("Processing precursors scans...");
		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));

		// get targeted ranges
		TDoubleHashSet boundaries=new TDoubleHashSet();
		ArrayList<Range> ranges=new ArrayList<Range>();
		for (Range range : stripefile.getRanges().keySet()) {
			boundaries.add(range.getStart());
			boundaries.add(range.getStop());
			ranges.add(range);
		}
		Collections.sort(ranges);
		
		double[] binBoundaries=boundaries.toArray();
		boolean[] useBin=new boolean[binBoundaries.length];
		Arrays.sort(binBoundaries);

		for (Range range : ranges) {
			int index=Arrays.binarySearch(binBoundaries, range.getMiddle());
			index=(-(index+1))-1;
			useBin[index]=true;
		}
		
		BlockingQueue<PeptideScoringResult> resultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
		PeptideScoringResultsConsumer resultsConsumer=taskFactory.getResultsConsumer(resultsQueue);
		Thread consumerThread=new Thread(resultsConsumer);
		consumerThread.start();
		
		int rangesFinished=0;
		// get stripes
		float numberOfTasks=2.0f+ranges.size();
		for (Range range : ranges) {
			String baseMessage="Working on "+range+" m/z";
			float baseIncrement=1.0f/numberOfTasks;
			float baseProgress=(1.0f+rangesFinished)/numberOfTasks;
			progress.update(baseMessage, baseProgress);
			int index=Arrays.binarySearch(binBoundaries, range.getMiddle());
			index=(-(index+1))-1;
			
			float dutyCycle=stripefile.getRanges().get(range);
			Logger.logLine("Processing "+range+" ("+dutyCycle+")");
			
			ArrayList<Stripe> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			Collections.sort(stripes);

			// prepare executor for background
			ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("SWATH_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
			LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
			ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

			int count=0;
			for (LibraryEntry entry : library.getEntries(range)) {
				count++;
				ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
				tasks.add(entry);
				tasks.add(entry.getReverse(parameters.getFragmentTolerance(), parameters.getAAConstants()));
				executor.submit(taskFactory.getScoringTask(pecanScorer, tasks, stripes, precursors, resultsQueue));
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

		consumerThread.join();
		resultsConsumer.close();

		progress.update("Running Percolator", (1.0f+rangesFinished)/numberOfTasks);
		ArrayList<ScoredObject<String>> passingPeptides=PercolatorExecutor.executePercolator(parameters.getPercolatorLocation(), featureFile, outputFile, parameters.getPercolatorThreshold());
		
		Logger.logLine("Finished analysis! "+resultsConsumer.getNumberProcessed()+" total peaks processed, "+passingPeptides.size()+" peaks identified at 1% FDR ("+(Math.round((System.currentTimeMillis()-startTime)/1000f/6f)/10f)+" minutes)");
		Logger.logLine(""); 
		progress.update(passingPeptides.size()+" peptides identified at "+(parameters.getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
	}

	public static boolean arePeptidesInRange(HashSet<FastaEntry> targets, Range range, PecanSearchParameters parameters) {
		// first check to see if we need to process this stripe
		boolean hasPeptides=false;
		outer:for (FastaEntry peptide : targets) {
			for (byte charge=parameters.getMinCharge(); charge<=parameters.getMaxCharge(); charge++) {
				double mz=parameters.getAAConstants().getChargedMass(peptide.getSequence(), charge);
				if (range.contains((float)mz)) {
					hasPeptides=true;
					break outer;
				}
			}
		}
		return hasPeptides;
	}

}
