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
import java.util.StringTokenizer;
import java.util.Map.Entry;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;

import com.google.common.base.Optional;
import com.google.common.util.concurrent.ThreadFactoryBuilder;

import edu.washington.gs.maccoss.encyclopedia.algorithms.BackgroundGenerator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.AbstractPecanFragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PecanLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.procedure.TDoubleObjectProcedure;
import gnu.trove.set.hash.TDoubleHashSet;

public class Pecanpie {
	public static void main(String[] args) {
		// EXAMPLE
		File diaFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/20150708_Ecoli_0911_25x4mzDIA_500_600.dia");
		//File fastaFile=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/ecoli-190209-contam_correctNL.fasta");
		File fastaFile=new File("/Users/searleb/Documents/projects/pecan/v0.9.7/ecoli_20150911_uniprot_sp_digested_Mass600to4000.fasta"); //FIXME
		File featureFile=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/encyc_report.feature.txt");
		File outputFile=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/encyc_report.percolator.txt");
		SearchParameters parameters=new SearchParameters(new AminoAcidConstants(), FragmentationType.YONLY, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"));
		PecanScoringFactory factory=new PecanOneScoringFactory(parameters, featureFile);
		
		ArrayList<FastaEntry> targets=new ArrayList<FastaEntry>();
		targets.add(new FastaEntry("FILE", ">Protein", "IGHTVEREDTPAIR"));
		//targets=null;
		
		try {
			runPie(Optional.fromNullable(targets), diaFile, fastaFile, featureFile, outputFile, factory);
		} catch (Exception e) {
			System.err.println("Encountered Fatal Error!");
			e.printStackTrace();
		}
	}

	public static void runPie(Optional<ArrayList<FastaEntry>> targetList, File diaFile, File fastaFile, File featureFile, File outputFile, PecanScoringFactory taskFactory) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		long startTime=System.currentTimeMillis();
		PSMScorer backgroundScorer=taskFactory.getBackgroundScorer();
		PSMScorer pecanScorer=taskFactory.getPecanScorer();
		SearchParameters parameters=taskFactory.getParameters();
		
		int cores=Runtime.getRuntime().availableProcessors();
		
		StripeFile stripefile=MzmlToDIAConverter.getFile(diaFile);

		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));

		HashSet<FastaEntry> targets=new HashSet<FastaEntry>();
		HashSet<String> backgroundProteome=new HashSet<String>();

		// pecan generates backgrounds using unique fasta peptides, target/decoy sequences, and 2000 random decoys for each window
		// add targets to proteome
		if (targetList.isPresent()) {
			for (FastaEntry target : targetList.get()) {
				targets.add(target);
				backgroundProteome.add(target.getSequence());
			}
		}
		
		// add database to proteome
		ArrayList<FastaEntry> entries=FastaReader.readFasta(fastaFile);
		for (FastaEntry entry : entries) {
			ArrayList<String> peptides=parameters.getEnzyme().digestProtein(entry.getSequence(), parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages());
			backgroundProteome.addAll(peptides);

			if (!targetList.isPresent()) {
				// search all peptides in database
				for (String peptide : peptides) {
					FastaEntry pe=entry.getSubEntry(peptide);
					targets.add(pe);
				}
			}
		}
		
		// get targeted ranges
		TDoubleHashSet boundaries=new TDoubleHashSet();
		ArrayList<Range> ranges=new ArrayList<Range>();
		for (Range range : stripefile.getRanges().keySet()) {
			boundaries.add(range.getStart());
			boundaries.add(range.getStop());
			if (arePeptidesInRange(targets, range, parameters)) {
				ranges.add(range);
			}
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
		
		Triplet<TDoubleIntHashMap[], ArrayList<String>[], HashSet<String>[]> background=BackgroundGenerator.generateBackground(binBoundaries, useBin, targets, backgroundProteome, parameters);
		TDoubleIntHashMap[] binCounters=background.x;
		ArrayList<String>[] backgroundProteomes=background.y;
		HashSet<String>[] backgroundDecoys=background.z;
		
		BlockingQueue<PeptideScoringResult> resultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
		PeptideScoringResultsConsumer resultsConsumer=taskFactory.getResultsConsumer(resultsQueue);
		Thread consumerThread=new Thread(resultsConsumer);
		consumerThread.start();
		
		// get stripes
		for (Range range : ranges) {
			int index=Arrays.binarySearch(binBoundaries, range.getMiddle());
			index=(-(index+1))-1;
			TDoubleIntHashMap map=binCounters[index];
			double[] keys=map.keys();
			Arrays.sort(keys);
			ArrayList<String> backgroundProteomeArray=backgroundProteomes[index];
			HashSet<String> backgroundProteomeSet=new HashSet<String>(backgroundProteomeArray);
			
			float dutyCycle=stripefile.getRanges().get(range);
			int scanAveragingMargin=(int)((parameters.getMinEluteTime())/dutyCycle); // floor
			float maxFragmentationMz=(float)Math.ceil(range.getMiddle()/10.0f)*20.0f+50.0f;
			Range fragmentationRange=new Range(maxFragmentationMz/15f, maxFragmentationMz);
			
			Logger.logLine("Processing "+range+" ("+scanAveragingMargin+")");
			
			ArrayList<Stripe> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			Collections.sort(stripes);

			// prepare executor for background
			ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("SWATH_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
			LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
			ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

			ArrayList<Future<HashMap<LibraryEntry, PeptideScoringResult>>> results=new ArrayList<Future<HashMap<LibraryEntry, PeptideScoringResult>>>();

			for (String peptide : backgroundDecoys[index]) {
				for (byte charge=parameters.getMinCharge(); charge<=parameters.getMaxCharge(); charge++) {
					double mz=parameters.getAAConstants().getChargedMass(peptide, charge);

					if (range.contains((float)mz)) {
						String random=PeptideUtils.getDecoy(peptide, backgroundProteomeSet, parameters);
						AbstractPecanFragmentationModel randmodel=taskFactory.getFragmentationModel(random, parameters.getAAConstants());
						PecanLibraryEntry randentry=randmodel.getPecanSpectrum(charge, keys, map, fragmentationRange, parameters, true);

						ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
						tasks.add(randentry);

						Future<HashMap<LibraryEntry, PeptideScoringResult>> value=executor.submit(new PeptideScoringTask(backgroundScorer, tasks, stripes, precursors));
						results.add(value);
					}
				}
			}
			executor.shutdown();
			while (!executor.isTerminated()) {
				Logger.logLine(workQueue.size()+" background peptides remaining for "+range+"...");
				Thread.sleep(500);
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
			results.clear();
			
			final ArrayList<XYPoint> means=new ArrayList<XYPoint>();
			final TDoubleObjectHashMap<XYPoint> backgroundScores=new TDoubleObjectHashMap<XYPoint>();
			backgroundScoreMap.forEachEntry(new TDoubleObjectProcedure<TDoubleArrayList>() {
				public boolean execute(double arg0, TDoubleArrayList arg1) {
					double[] values=arg1.toArray();
					double m=General.mean(values);
					double s=General.stdev(values);
					backgroundScores.put(arg0, new XYPoint(m, s));
					means.add(new XYPoint(arg0, m));
					return true;
				};
			});
			//Charter.launchChart("RT ("+range+" M/Z)", "Fragment Intensity", true, new XYTrace(means, GraphType.line, "Background"));

			for (FastaEntry peptide : targets) {
				String sequence=peptide.getSequence();
				for (byte charge=parameters.getMinCharge(); charge<=parameters.getMaxCharge(); charge++) {
					double mz=parameters.getAAConstants().getChargedMass(sequence, charge);
					if (range.contains((float)mz)) {
						AbstractPecanFragmentationModel model=taskFactory.getFragmentationModel(sequence, parameters.getAAConstants());
						PecanLibraryEntry pecanEntry=model.getPecanSpectrum(charge, keys, map, fragmentationRange, parameters, false);
						
						AbstractPecanFragmentationModel revmodel=taskFactory.getFragmentationModel(PeptideUtils.getSmartDecoy(sequence, charge, backgroundProteomeSet, parameters), parameters.getAAConstants());
						PecanLibraryEntry reventry=revmodel.getPecanSpectrum(charge, keys, map, fragmentationRange, parameters, true);

						ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
						tasks.add(pecanEntry);
						tasks.add(reventry);

						executor.submit(taskFactory.getScoringTask(pecanScorer, tasks, stripes, backgroundScores, precursors, scanAveragingMargin, resultsQueue));
					}
				}
			}
			executor.shutdown();
			while (!executor.isTerminated()) {
				Logger.logLine(workQueue.size()+" peptides remaining for "+range+"...");
				Thread.sleep(500);
			}
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
		}
		resultsQueue.put(PeptideScoringResult.POISON_RESULT);

		consumerThread.join();
		resultsConsumer.close();
		
		PercolatorExecutor e=new PercolatorExecutor(featureFile);
		BlockingQueue<OutputMessage> result=e.start();
		
		int countBelow1pFDR=0;
		boolean isFirst=true;
		boolean record=true;
		PrintWriter writer=new PrintWriter(outputFile, "UTF-8");
		while (!e.isFinished()||!result.isEmpty()) {
			if (!result.isEmpty()) {
				OutputMessage data=result.take();
				if (data.isStdOutput) {
					if (isFirst) {
						isFirst=false;
					} else if (record) {
						StringTokenizer st=new StringTokenizer(data.message);
						st.nextToken(); // PSMid
						st.nextToken(); // score
						float qvalue=Float.parseFloat(st.nextToken());
						if (qvalue<0.01f) {
							countBelow1pFDR++;
						} else {
							record=false;
						}
					}
					writer.println(data.message);
				} else {
					Logger.logLine(data.message);
				}
			} else {
				Thread.sleep(10);
			}
		}
		writer.flush();
		writer.close();
		
		Logger.logLine("Finished analysis! "+resultsConsumer.getNumberProcessed()+" total peaks processed, "+countBelow1pFDR+" peaks identified at 1% FDR ("+(Math.round((System.currentTimeMillis()-startTime)/1000f/6f)/10f)+" minutes)");
	}

	public static boolean arePeptidesInRange(HashSet<FastaEntry> targets, Range range, SearchParameters parameters) {
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
