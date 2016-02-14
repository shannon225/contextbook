package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

import edu.washington.gs.maccoss.encyclopedia.algorithms.DotProduct;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringResultsToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchFeatureReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TDoubleHashSet;

public class ChromatogramExtractor {
	public static void parse(ArrayList<ScoredObject<String>> globalPassingPeptides, SearchJobData job, LibraryFile libraryFile) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		int minChromatogramCount=4;
		SearchParameters parameters=job.getParameters();
		PSMScorer scorer=new DotProduct(parameters.getFragmentTolerance());
		float duration=3*parameters.getExpectedPeakWidth(); // search for 6 minutes
		
		ArrayList<ScoredObject<String>>	localPassingPeptides=PercolatorReader.getPassingPeptidesFromTSV(job.getOutputFile(), parameters.getPercolatorThreshold());
		HashSet<String> passingPeptideSequences=new HashSet<String>();
		for (ScoredObject<String> psm : globalPassingPeptides) {
			String peptideModSeq=PecanScoringResultsToTSVConsumer.getPeptideSequence(psm.y);
			passingPeptideSequences.add(peptideModSeq);
		}
		
		ArrayList<LibraryEntry> entries=libraryFile.getAllEntries(false);
		HashMap<String, LibraryEntry> entireEntryMap=new HashMap<String, LibraryEntry>();
		for (LibraryEntry entry : entries) {
			LibraryEntry unitEntry=entries.get(0).toUnitSpectrum(minChromatogramCount);
			if (unitEntry.getIonCount()<minChromatogramCount) {
				entireEntryMap.put(entry.getPeptideModSeq(), unitEntry);
			}
		}
		
		final HashMap<String, LibraryEntry> keptEntryMap=new HashMap<String, LibraryEntry>();
		for (ScoredObject<String> psm : localPassingPeptides) {
			String peptideModSeq=PecanScoringResultsToTSVConsumer.getPeptideSequence(psm.y);
			LibraryEntry entry=entireEntryMap.get(peptideModSeq);
			if (entry!=null&&passingPeptideSequences.contains(peptideModSeq)) {
				keptEntryMap.put(peptideModSeq, entry);
			}
		}
				
		File diaFile=job.getDiaFile();
		StripeFileInterface stripefile=MzmlToDIAConverter.getFile(diaFile, job.getParameters());
		File featureFile=job.getFeatureFile();
		HashMap<String, PSM> psmData=getPSMData(featureFile);

		TDoubleHashSet boundaries=new TDoubleHashSet();
		ArrayList<Range> ranges=new ArrayList<Range>();
		for (Range range : stripefile.getRanges().keySet()) {
			boundaries.add(range.getStart());
			boundaries.add(range.getStop());
			if (!parameters.useTargetWindowCenter()||range.contains(parameters.getTargetWindowCenter())) {
				ranges.add(range);
			}
		}
		Collections.sort(ranges);

		File integrationFile=new File(diaFile.getAbsolutePath()+".quant.txt");
		PrintWriter writer=new PrintWriter(integrationFile, "UTF-8");
		writer.println("File\tPeptideModSeq\tPrecursorCharge\tFragmentIons\tRTStart\tRTCenter\tRTStop\tTIC\tNormTIC");

		for (Range range : ranges) {
			Logger.logLine("Working on "+range+" m/z");

			ArrayList<Stripe> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			TIntArrayList stripeScanIDs=new TIntArrayList();
			for (Stripe stripe : stripes) {
				stripeScanIDs.add(stripe.getSpectrumIndex());
			}
			Collections.sort(stripes);

			ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("SWATH_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
			LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
			ExecutorService executor=new ThreadPoolExecutor(parameters.getNumberOfThreadsUsed(), parameters.getNumberOfThreadsUsed(), Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

			for (LibraryEntry unitEntry : keptEntryMap.values()) {
				PSM psm=psmData.get(unitEntry.getPeptideModSeq());
				
				if (range.contains((float)psm.precursorMZ)) {
					int index=stripeScanIDs.binarySearch(psm.scanID);
					if (index<0) {
						if (index<0) {
							index=-(index+1);
						}
					}

					ArrayList<Stripe> group=new ArrayList<Stripe>();
					for (int i=index; i<stripes.size(); i++) {
						if (stripes.get(i).getScanStartTime()<psm.scanStartTime+duration) {
							group.add(stripes.get(i));
						}
					}
					for (int i=index-1; i>=0; i--) {
						if (stripes.get(i).getScanStartTime()>psm.scanStartTime-duration) {
							group.add(0, stripes.get(i));
						}
					}
					Triplet<double[], float[], Range> result=SearchFeatureReader.quantifyPeptide(scorer, unitEntry, psm.peptideModSeq, psm.scanStartTime, true, stripes);
					float[] intensities=result.y;
					float sum=0.0f;
					for (int i=0; i<intensities.length; i++) {
						sum+=intensities[i];
					}
					writer.println(diaFile.getName()+"\t"+psm.peptideModSeq+"\t"+psm.precursorCharge+"\t"+result.z.getStart()+"\t"+psm.scanStartTime+"\t"+result.z.getStop()+"\t"+sum);
				}
			}
			
			executor.shutdown();
			while (!executor.isTerminated()) {
				Logger.logLine(workQueue.size()+" peptides remaining for "+range+"...");
				Thread.sleep(500);
			}
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
		}
		
		writer.flush();
		writer.close();
	}
	
	private static HashMap<String, PSM> getPSMData(File featureFile) throws InterruptedException {
		final HashMap<String, PSM> map=new HashMap<String, ChromatogramExtractor.PSM>();
		BlockingQueue<Map<String, String>> blockingQueue=new LinkedBlockingQueue<Map<String, String>>();
		TableParserProducer producer=new TableParserProducer(blockingQueue, featureFile, "\t", 1);
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String psmID=row.get("id");
				int scanID=Integer.parseInt(row.get("ScanNr"));
				double precursorMZ=Double.parseDouble(row.get("precursorMz"));
				byte precursorCharge=PecanScoringResultsToTSVConsumer.getCharge(psmID);
				String peptideModSeq=PecanScoringResultsToTSVConsumer.getPeptideSequence(psmID);
				String rtString=row.get("midTime");
				if (rtString==null) rtString=row.get("RTinMin");
				float retentionTime=Float.parseFloat(rtString)*60f;
				map.put(psmID, new PSM(scanID, precursorMZ, precursorCharge, retentionTime, peptideModSeq));
			}
		};
		TableParserConsumer consumer=new TableParserConsumer(blockingQueue, muscle);

		Thread producerThread=new Thread(producer);
		producerThread.start();

		Thread consumerThread=new Thread(consumer);
		consumerThread.start();

		producerThread.join();
		consumerThread.join();
		return map;
	}
	
	private static class PSM {
		private final int scanID;
		private final double precursorMZ;
		private final byte precursorCharge;
		private final String peptideModSeq;
		private final float scanStartTime;
		public PSM(int scanID, double precursorMZ, byte precursorCharge, float scanStartTime, String peptideModSeq) {
			this.scanID=scanID;
			this.precursorMZ=precursorMZ;
			this.precursorCharge=precursorCharge;
			this.scanStartTime=scanStartTime;
			this.peptideModSeq=peptideModSeq;
		}
		
	}
}
