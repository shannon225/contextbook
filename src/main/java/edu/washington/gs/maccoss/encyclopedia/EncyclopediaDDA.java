package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackground;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.SaveResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.TeeResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;
import java.util.zip.DataFormatException;

public class EncyclopediaDDA {
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
		ArrayList<Range> ranges=library.getMinMaxMZ().chunkIntoBins(200);

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

		// prepare executor for background
		ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("LIBRARY-%d").setDaemon(true).build();
		LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
		ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

		int count=0;
		// get stripes
		for (Range range : ranges) {
			Logger.logLine("Processing "+range);
			
			List<Stripe> stripes=stripefile.getStripes(range, -Float.MAX_VALUE, Float.MAX_VALUE, true);
			ArrayList<LibraryEntry> entries=library.getEntries(range, true);
			LibraryBackgroundInterface background=new LibraryBackground(entries);
			PSMScorer scorer=taskFactory.getLibraryScorer(background);
			ArrayList<LibraryEntry> reverses=new ArrayList<LibraryEntry>();
			for (LibraryEntry entry : entries) {
				count++;
				reverses.add(entry.getDecoy(parameters));

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
					reverses.add(shuffle);
					reverses.add(shuffle.getDecoy(parameters));
				}
				
			}
			entries.addAll(reverses);
			executor.submit(taskFactory.getDDAScoringTask(scorer, entries, stripes, precursors, resultsQueue));
		}
		executor.shutdown();
		while (!executor.isTerminated()) {
			Logger.logLine(workQueue.size()+" / "+ranges.size()+" windows remaining");
			float finishedFraction=(ranges.size()-workQueue.size())/(float)ranges.size();
			progress.update(workQueue.size()+" / "+ranges.size()+" windows remaining", finishedFraction*0.9f);
			Thread.sleep(500);
		}
		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
		
		resultsQueue.put(PeptideScoringResult.POISON_RESULT);

		consumer1Thread.join();
		consumer2Thread.join();
		consumer3Thread.join();
		teeConsumer.close();
		progress.update("Organizing results", 0.9f);

		ArrayList<PercolatorPeptide> passingPeptides=Encyclopedia.percolatePeptides(progress, job, stripefile, saveResultsConsumer);
		
		Logger.logLine("Finished analysis! "+writeResultsConsumer.getNumberProcessed()+" total peaks processed, "+passingPeptides.size()+" peaks identified at "+(parameters.getEffectivePercolatorThreshold()*100f)+"% FDR ("+(Math.round((System.currentTimeMillis()-startTime)/1000f/6f)/10f)+" minutes)");
		Logger.logLine(""); 
	}
}
