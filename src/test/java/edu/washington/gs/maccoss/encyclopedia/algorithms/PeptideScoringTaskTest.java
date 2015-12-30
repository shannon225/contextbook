package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneFragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanRawScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PecanLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.set.hash.TDoubleHashSet;

public class PeptideScoringTaskTest {
	private static final byte PLOTTING_METHOD=FragmentationTraceTask.PLOT_INTENSITIES;
	
	private static final Range FRAGMENTATION_RANGE=new Range(0f, 2000f);
	private static final SearchParameters PARAMETERS=new SearchParameters(new AminoAcidConstants(), FragmentationType.YONLY, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"));

	public static void main(String[] args) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		int cores=Runtime.getRuntime().availableProcessors();

		File f=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/20150708_Ecoli_0931_25x4mzDIA_700_800.dia");
		StripeFileInterface stripefile=new StripeFile();
		stripefile.openFile(f);

		byte[] charges=new byte[] {(byte)3};
		String[] peptides=new String[] {"QHLYSISDEQLRPYFPENK"};

		InputStream is=stripefile.getClass().getResourceAsStream("/ecoli-190209-contam_correctNL.fasta");
		ArrayList<FastaEntry> entries=FastaReader.readFasta(is, "ecoli-190209-contam_correctNL.fasta");
		
		TDoubleHashSet boundaries=new TDoubleHashSet();
		for (Range range : stripefile.getRanges().keySet()) {
			boundaries.add(range.getStart());
			boundaries.add(range.getStop());
		}
		double[] binArray=boundaries.toArray();
		Arrays.sort(binArray);
		Pair<TDoubleIntHashMap[], ArrayList<String>[]> background=BackgroundGenerator.generateBackground(binArray, entries, true, PARAMETERS);
		TDoubleIntHashMap[] binCounters=background.x;
		ArrayList<String>[] backgroundProteomes=background.y;
		
		PecanRawScorer pecanScorer=new PecanRawScorer(PARAMETERS.getFragmentTolerance(), new ExpectedFragmentationScorer(PARAMETERS));
		
		// get stripes
		for (Entry<Range, Float> entry : stripefile.getRanges().entrySet()) {
			Range range=entry.getKey();
			float dutyCycle=entry.getValue();
			int scanAveragingMargin=(int)(PARAMETERS.getMinEluteTime()/dutyCycle/2); // floor
			
			System.out.println("Processing "+range+" ("+scanAveragingMargin+")");
			
			ArrayList<Stripe> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			Collections.sort(stripes);
			
			int index=Arrays.binarySearch(binArray, range.getMiddle());
			index=(-(index+1))-1;
			TDoubleIntHashMap map=binCounters[index];
			double[] keys=map.keys();
			Arrays.sort(keys);
			ArrayList<String> backgroundProteomeArray=backgroundProteomes[index];
			
			if (peptides!=null) {
				backgroundProteomeArray=new ArrayList<String>(Arrays.asList(peptides));
			}
			
			// first check to see if we need to process this stripe
			boolean hasPeptides=false;
			outer:for (String peptide : backgroundProteomeArray) {
				for (byte charge : charges) {
					double mz=PARAMETERS.getAAConstants().getChargedMass(peptide, charge);
					if (range.contains((float)mz)) {
						hasPeptides=true;
						break outer;
					}
				}
			}
			if (!hasPeptides) {
				continue;
			}

			// prepare executor for background
			ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("SWATH_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
			LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
			ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

			ArrayList<Future<HashMap<LibraryEntry, PeptideScoringResult>>> results=new ArrayList<Future<HashMap<LibraryEntry, PeptideScoringResult>>>();
			for (String peptide : backgroundProteomeArray) {
				for (byte charge : charges) {
					double mz=PARAMETERS.getAAConstants().getChargedMass(peptide, charge);
					if (range.contains((float)mz)) {
						PecanOneFragmentationModel model=new PecanOneFragmentationModel(peptide, PARAMETERS.getAAConstants());
						PecanLibraryEntry pecanEntry=model.getPecanSpectrum(charge, keys, map, FRAGMENTATION_RANGE, PARAMETERS, false);

						ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
						tasks.add(pecanEntry);

						Future<HashMap<LibraryEntry, PeptideScoringResult>> value=executor.submit(new FragmentationTraceTask(pecanScorer, PLOTTING_METHOD, tasks, stripes, new PrecursorScanMap(new ArrayList<PrecursorScan>())));
						results.add(value);
					}
				}
			}
			executor.shutdown();
			while (!executor.isTerminated()) {
				System.out.println(workQueue.size()+" peptides remaining for "+range+"...");
				Thread.sleep(100);
			}
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

			ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
			for (Future<HashMap<LibraryEntry, PeptideScoringResult>> future : results) {
				HashMap<LibraryEntry, PeptideScoringResult> result=future.get();
				for (Entry<LibraryEntry, PeptideScoringResult> resultEntry : result.entrySet()) {
					//LibraryEntry peptide=resultEntry.getKey();
					FragmentationScoringResult peptideResult=(FragmentationScoringResult)resultEntry.getValue();

					for (XYTrace trace : peptideResult.getFragmentationTraces()) {
						traces.add(trace);
					}
				}
			}
			String yAxisName;
			if (FragmentationTraceTask.PLOT_INTENSITIES==PLOTTING_METHOD) {
				yAxisName="Intensity";
			} else if (FragmentationTraceTask.PLOT_SCORES==PLOTTING_METHOD) {
				yAxisName="Score";
			} else if (FragmentationTraceTask.PLOT_DELTA_MASSES==PLOTTING_METHOD) {
				yAxisName="Delta Mass";

			}
			Charter.launchChart("RT ("+range+" M/Z)", yAxisName, true, traces.toArray(new XYTrace[traces.size()]));
		}
	}


}
