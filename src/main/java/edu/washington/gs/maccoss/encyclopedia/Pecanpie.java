package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
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

import edu.washington.gs.maccoss.encyclopedia.algorithms.BackgroundGenerator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.DotProduct;
import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanAuxillaryScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanRawScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringTask;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PecanLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.procedure.TDoubleObjectProcedure;
import gnu.trove.set.hash.TDoubleHashSet;

public class Pecanpie {
	public static void main(String[] args) {
		// EXAMPLE
		File diaFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/20150708_Ecoli_0911_25x4mzDIA_500_600.dia");
		File fastaFile=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/ecoli-190209-contam_correctNL.fasta");
		File outputFile=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/encyc_report.txt");
		SearchParameters parameters=new SearchParameters(new AminoAcidConstants(), FragmentationType.YONLY, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"));

		PSMScorer backgroundScorer=new DotProduct(parameters.getFragmentTolerance());
		PSMScorer pecanScorer=new PecanRawScorer(parameters.getFragmentTolerance(), new PecanAuxillaryScorer(parameters));
		
		try {
			runPie(diaFile, fastaFile, outputFile, backgroundScorer, pecanScorer, parameters);
		} catch (Exception e) {
			System.err.println("Encountered Fatal Error!");
			e.printStackTrace();
		}
	}

	public static void runPie(File diaFile, File fastaFile, File outputFile, PSMScorer backgroundScorer, PSMScorer pecanScorer, SearchParameters parameters) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		int cores=Runtime.getRuntime().availableProcessors();
		
		StripeFile stripefile=new StripeFile();
		stripefile.openFile(diaFile);

		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));

		ArrayList<FastaEntry> entries=FastaReader.readFasta(fastaFile);

		PrintWriter writer = new PrintWriter(outputFile, "UTF-8");
		
		TDoubleHashSet boundaries=new TDoubleHashSet();
		ArrayList<Range> ranges=new ArrayList<Range>();
		for (Range range : stripefile.getRanges().keySet()) {
			boundaries.add(range.getStart());
			boundaries.add(range.getStop());
			ranges.add(range);
		}
		Collections.sort(ranges);
		
		double[] binArray=boundaries.toArray();
		Arrays.sort(binArray);
		
		Pair<TDoubleIntHashMap[], ArrayList<String>[]> background=BackgroundGenerator.generateBackground(binArray, entries, parameters);
		TDoubleIntHashMap[] binCounters=background.x;
		ArrayList<String>[] backgroundProteomes=background.y;

		byte[] charges=new byte[parameters.getMaxCharge()-parameters.getMinCharge()+1];
		for (int i=0; i<charges.length; i++) {
			charges[i]=(byte)(parameters.getMinCharge()+i);
		}
		
		// get stripes
		for (Range range : ranges) {
			float dutyCycle=stripefile.getRanges().get(range);
			int scanAveragingMargin=(int)(((parameters.getMinEluteTime())/dutyCycle+1)/2); // floor
			
			System.out.println("Processing "+range+" ("+scanAveragingMargin+")");
			
			int index=Arrays.binarySearch(binArray, range.getMiddle());
			index=(-(index+1))-1;
			TDoubleIntHashMap map=binCounters[index];
			double[] keys=map.keys();
			Arrays.sort(keys);
			ArrayList<String> backgroundProteomeArray=backgroundProteomes[index];
			HashSet<String> backgroundProteomeSet=new HashSet<String>(backgroundProteomeArray);
			
			ArrayList<String> targetPeptides=backgroundProteomeArray;
			
			// first check to see if we need to process this stripe
			boolean hasPeptides=false;
			outer:for (String peptide : targetPeptides) {
				for (byte charge : charges) {
					double mz=parameters.getAAConstants().getChargedMass(peptide, charge);
					if (range.contains((float)mz)) {
						hasPeptides=true;
						break outer;
					}
				}
			}
			if (!hasPeptides) {
				continue;
			}
			
			ArrayList<Stripe> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			Collections.sort(stripes);

			// prepare executor for background
			ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("SWATH_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
			LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
			ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

			int backgroundPeptideCount=0;
			int seed=RandomGenerator.randomInt(1);
			ArrayList<Future<HashMap<LibraryEntry, PeptideScoringResult>>> results=new ArrayList<Future<HashMap<LibraryEntry, PeptideScoringResult>>>();
			while (backgroundPeptideCount<2000) {
				for (byte charge : charges) {
					seed=RandomGenerator.randomInt(seed);
					String peptide=backgroundProteomeArray.get((int)(RandomGenerator.floatFromRandomInt(seed)*backgroundProteomeArray.size()));
					double mz=parameters.getAAConstants().getChargedMass(peptide, charge);

					if (range.contains((float)mz)) {
						String random=PeptideUtils.getDecoy(peptide, backgroundProteomeSet, parameters);
						FragmentationModel randmodel=new FragmentationModel(random, parameters.getAAConstants());
						PecanLibraryEntry randentry=randmodel.getPecanSpectrum(charge, keys, map, parameters, true);

						ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
						tasks.add(randentry);

						Future<HashMap<LibraryEntry, PeptideScoringResult>> value=executor.submit(new PeptideScoringTask(backgroundScorer, tasks, stripes, precursors));
						results.add(value);
						
						backgroundPeptideCount++;
					}
				}
			}
			executor.shutdown();
			while (!executor.isTerminated()) {
				System.out.println(workQueue.size()+" background peptides remaining for "+range+"...");
				Thread.sleep(200);
			}
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

			// prepare executor for peptides
			executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 
			
			TDoubleObjectHashMap<TDoubleArrayList> backgroundScoreMap=new TDoubleObjectHashMap<TDoubleArrayList>();
			for (Future<HashMap<LibraryEntry, PeptideScoringResult>> future : results) {
				HashMap<LibraryEntry, PeptideScoringResult> result=future.get();
				for (Entry<LibraryEntry, PeptideScoringResult> resultEntry : result.entrySet()) {
					Pair<double[], double[]> arrays=resultEntry.getValue().getTrace().toArrays();
					double[] x=arrays.x;
					double[] y=arrays.y;
					for (int i=0; i<x.length; i++) {
						TDoubleArrayList list=backgroundScoreMap.get(x[i]);
						if (list==null) {
							list=new TDoubleArrayList();
							backgroundScoreMap.put(x[i], list);
						}
						list.add(y[i]);
					}
				}
			}
			
			final TDoubleObjectHashMap<XYPoint> backgroundScores=new TDoubleObjectHashMap<XYPoint>();
			backgroundScoreMap.forEachEntry(new TDoubleObjectProcedure<TDoubleArrayList>() {
				public boolean execute(double arg0, TDoubleArrayList arg1) {
					double[] values=arg1.toArray();
					double m=General.mean(values);
					double s=General.stdev(values);
					backgroundScores.put(arg0, new XYPoint(m, s));
					return true;
				};
			});

			results.clear();
			for (String peptide : targetPeptides) {
				for (byte charge : charges) {
					double mz=parameters.getAAConstants().getChargedMass(peptide, charge);
					if (range.contains((float)mz)) {
						FragmentationModel model=new FragmentationModel(peptide, parameters.getAAConstants());
						PecanLibraryEntry pecanEntry=model.getPecanSpectrum(charge, keys, map, parameters, false);

						FragmentationModel revmodel=new FragmentationModel(PeptideUtils.getSmartDecoy(peptide, charge, backgroundProteomeSet, parameters), parameters.getAAConstants());
						PecanLibraryEntry reventry=revmodel.getPecanSpectrum(charge, keys, map, parameters, true);

						ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
						tasks.add(pecanEntry);
						tasks.add(reventry);

						Future<HashMap<LibraryEntry, PeptideScoringResult>> value=executor.submit(new PecanScoringTask(pecanScorer, tasks, stripes, backgroundScores, precursors, scanAveragingMargin));
						results.add(value);
					}
				}
			}
			executor.shutdown();
			while (!executor.isTerminated()) {
				System.out.println(workQueue.size()+" peptides remaining for "+range+"...");
				Thread.sleep(200);
			}
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

			ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
			for (Future<HashMap<LibraryEntry, PeptideScoringResult>> future : results) {
				HashMap<LibraryEntry, PeptideScoringResult> result=future.get();
				for (Entry<LibraryEntry, PeptideScoringResult> resultEntry : result.entrySet()) {
					PecanLibraryEntry peptide=(PecanLibraryEntry)resultEntry.getKey();
					PeptideScoringResult peptideResult=resultEntry.getValue();
					
					int rank=1;
					for (Pair<ScoredObject<Stripe>, float[]> goodStripe : peptideResult.getGoodStripes()) {
						float primaryScore=goodStripe.x.x;
						Stripe stripe=goodStripe.x.y;
						float[] auxScores=goodStripe.y;
						
						if (rank<=3) {
							writer.print(peptide.getPeptideModSeq()+"\t"+peptide.isDecoy()+"\t"+rank+"\t"+primaryScore+"\t"+stripe.getScanStartTime());
							for (float s : auxScores) {
								writer.print("\t"+s);
							}
							writer.println();
						}
						rank++;
						if (rank>3) break;
					}
					traces.add(peptideResult.getTrace());
				}
			}
			writer.flush();
		}
		writer.close();
	}

}
