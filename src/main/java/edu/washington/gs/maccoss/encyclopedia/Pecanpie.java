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
import java.util.StringTokenizer;
import java.util.TreeMap;
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
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.AbstractPecanFragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PecanLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.procedure.TDoubleObjectProcedure;
import gnu.trove.set.hash.TDoubleHashSet;

public class Pecanpie {
	public static void main(String[] args) {
		Logger.logLine("Pecanpie version 0.1");
		
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		if (arguments.size()==0||arguments.containsKey("-v")||arguments.containsKey("-version")||arguments.containsKey("--version")||arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.logLine("Required Parameters: ");
			Logger.logLine("\t-i\tinput .DIA or .MZML file");
			Logger.logLine("\t-f\tbackground FASTA file");
			Logger.logLine("Other Parameters: ");
			Logger.logLine("\t-t\ttarget FASTA file (default: background FASTA file)");
			Logger.logLine("\t-o\toutput report file (default: [input file].pecan.txt)");
			
			TreeMap<String, String> defaults=new TreeMap<String, String>(SearchParameterParser.getDefaultParameters());
			for (Entry<String, String> entry : defaults.entrySet()) {
				Logger.logLine("\t"+entry.getKey()+"\t(default: "+entry.getValue()+")");
			}
			System.exit(1);
		}
		
		if (!arguments.containsKey("-i")||!arguments.containsKey("-f")) {
			Logger.errorLine("You are required to specify an input file (-i) and a background FASTA file (-f)");
			System.exit(1);
		}

		File diaFile=new File(arguments.get("-i"));
		File fastaFile=new File(arguments.get("-f"));
		
		ArrayList<FastaEntry> targets;
		if (arguments.containsKey("-t")) {
			targets=FastaReader.readFasta(new File(arguments.get("-t")));
		} else {
			targets=null;
		}
		
		File outputFile;
		if (arguments.containsKey("-o")) {
			outputFile=new File(arguments.get("-o"));
		} else {
			outputFile=new File(diaFile.getAbsolutePath()+".pecan.txt");
		}
		
		
		File featureFile=new File(outputFile.getAbsolutePath()+".features.txt");

		SearchParameters parameters=SearchParameterParser.parseParameters(arguments);
		PecanScoringFactory factory=new PecanOneScoringFactory(parameters, featureFile);

		Logger.logLine("Parameters:");
		Logger.logLine(" -i "+diaFile.getAbsolutePath());
		Logger.logLine(" -f "+fastaFile.getAbsolutePath());
		Logger.logLine(" -t "+arguments.get("-t"));
		Logger.logLine(" -o "+outputFile.getAbsolutePath());
		Logger.logLine(parameters.toString());
		
		try {
			runPie(new EmptyProgressIndicator(), Optional.fromNullable(targets), diaFile, fastaFile, featureFile, outputFile, factory);
		} catch (Exception e) {
			System.err.println("Encountered Fatal Error!");
			e.printStackTrace();
		}
	}

	public static void runPie(ProgressIndicator progress, Optional<ArrayList<FastaEntry>> targetList, File diaFile, File fastaFile, File featureFile, File outputFile, PecanScoringFactory taskFactory) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		long startTime=System.currentTimeMillis();
		PSMScorer backgroundScorer=taskFactory.getBackgroundScorer();
		PSMScorer pecanScorer=taskFactory.getPecanScorer();
		SearchParameters parameters=taskFactory.getParameters();
		
		int cores=Runtime.getRuntime().availableProcessors();

		progress.update("Converting files...", Float.MIN_VALUE);
		
		StripeFileInterface stripefile=MzmlToDIAConverter.getFile(diaFile);

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
		
		int rangesFinished=0;
		// get stripes
		for (Range range : ranges) {
			progress.update("Working on "+range+" m/z", (1.0f+rangesFinished)/(2.0f+ranges.size()));
			int index=Arrays.binarySearch(binBoundaries, range.getMiddle());
			index=(-(index+1))-1;
			TDoubleIntHashMap map=binCounters[index];
			double[] keys=map.keys();
			Arrays.sort(keys);
			ArrayList<String> backgroundProteomeArray=backgroundProteomes[index];
			HashSet<String> backgroundProteomeSet=new HashSet<String>(backgroundProteomeArray);
			
			float dutyCycle=stripefile.getRanges().get(range);
			int scanAveragingMargin=(int)((parameters.getMinEluteTime())/dutyCycle)+1; // floor
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
						ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
						
						AbstractPecanFragmentationModel model=taskFactory.getFragmentationModel(sequence, parameters.getAAConstants());
						PecanLibraryEntry pecanEntry=model.getPecanSpectrum(charge, keys, map, fragmentationRange, parameters, false);
						tasks.add(pecanEntry);
						
						if (!parameters.isDontRunDecoys()) {
							AbstractPecanFragmentationModel revmodel=taskFactory.getFragmentationModel(PeptideUtils.getSmartDecoy(sequence, charge, backgroundProteomeSet, parameters), parameters.getAAConstants());
							PecanLibraryEntry reventry=revmodel.getPecanSpectrum(charge, keys, map, fragmentationRange, parameters, true);
							tasks.add(reventry);
						}

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
			
			rangesFinished++;
		}
		resultsQueue.put(PeptideScoringResult.POISON_RESULT);

		consumerThread.join();
		resultsConsumer.close();

		progress.update("Running Percolator", (1.0f+rangesFinished)/(2.0f+ranges.size()));
		PercolatorExecutor e=new PercolatorExecutor(featureFile);
		BlockingQueue<OutputMessage> result=e.start();
		
		int countBelow1pFDR=0;
		boolean isFirst=true;
		boolean record=true;
		PrintWriter writer=new PrintWriter(outputFile, "UTF-8");
		while (!e.isFinished()||!result.isEmpty()) {
			if (!result.isEmpty()) {
				OutputMessage data=result.take();
				if (data.isStdOutput()) {
					if (isFirst) {
						isFirst=false;
					} else if (record) {
						StringTokenizer st=new StringTokenizer(data.getMessage());
						st.nextToken(); // PSMid
						st.nextToken(); // score
						float qvalue=Float.parseFloat(st.nextToken());
						if (qvalue<0.01f) {
							countBelow1pFDR++;
						} else {
							record=false;
						}
					}
					writer.println(data.getMessage());
				} else {
					Logger.logLine(data.getMessage());
				}
			} else {
				Thread.sleep(10);
			}
		}
		writer.flush();
		writer.close();
		
		Logger.logLine("Finished analysis! "+resultsConsumer.getNumberProcessed()+" total peaks processed, "+countBelow1pFDR+" peaks identified at 1% FDR ("+(Math.round((System.currentTimeMillis()-startTime)/1000f/6f)/10f)+" minutes)");
		progress.update(countBelow1pFDR+" peaks identified at 1% FDR", 1.0f);
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
