package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import edu.washington.gs.maccoss.encyclopedia.algorithms.*;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.*;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
import edu.washington.gs.maccoss.encyclopedia.filereaders.*;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.*;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.procedure.TDoubleObjectProcedure;
import gnu.trove.set.hash.TDoubleHashSet;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.*;
import java.util.Map.Entry;
import java.util.concurrent.*;
import java.util.function.Supplier;
import java.util.zip.DataFormatException;

public class Pecanpie {
	public static final String TARGET_FASTA_TAG="-t";
	public static final String OUTPUT_RESULT_TAG="-o";
	public static final String INPUT_DIA_TAG="-i";
	public static final String BACKGROUND_FASTA_TAG="-f";

	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		if (arguments.size()==0) {
			SearchGUIMain.runGUI(ProgramType.PecanPie);
			
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("Walnut Help");
			Logger.timelessLogLine("Walnut is a FASTA database search engine for DIA data that uses PECAN-style scoring.");
			Logger.timelessLogLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .DIA or .MZML file");
			Logger.timelessLogLine("\t-f\tbackground FASTA file");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-t\ttarget FASTA file (default: background FASTA file)");
			Logger.timelessLogLine("\t-tp\ttrue/false target FASTA file contains peptides (default: false)"); 
			Logger.timelessLogLine("\t-o\toutput report file (default: [input file]"+PecanJobData.OUTPUT_FILE_SUFFIX+")");
			
			TreeMap<String, String> defaults=new TreeMap<String, String>(PecanParameterParser.getDefaultParameters());
			int maxWidth=0;
			for (String key : defaults.keySet()) {
				if (key.length()>maxWidth) maxWidth=key.length();
			}
			for (Entry<String, String> entry : defaults.entrySet()) {
				Logger.timelessLogLine("\t"+General.formatCellToWidth(entry.getKey(), maxWidth)+" (default: "+entry.getValue()+")");
			}
			System.exit(1);
			
		} else if (arguments.containsKey("-v")||arguments.containsKey("-version")||arguments.containsKey("--version")) {
			Logger.logLine("Walnut version "+PecanOneScoringFactory.version);
			System.exit(1);
			
		} else {
			if (!arguments.containsKey(INPUT_DIA_TAG)||!arguments.containsKey(BACKGROUND_FASTA_TAG)) {
				Logger.errorLine("You are required to specify an input file ("+INPUT_DIA_TAG+") and a background FASTA file ("+BACKGROUND_FASTA_TAG+")");
				System.exit(1);
			}

			File diaFile=new File(arguments.get(INPUT_DIA_TAG));
			File fastaFile=new File(arguments.get(BACKGROUND_FASTA_TAG));

			File outputFile;
			if (arguments.containsKey(OUTPUT_RESULT_TAG)) {
				outputFile=new File(arguments.get(OUTPUT_RESULT_TAG));
			} else {
				outputFile=new File(diaFile.getAbsolutePath()+PecanJobData.OUTPUT_FILE_SUFFIX);
			}

			File decoyFile=new File(PecanJobData.getOutputAbsolutePathPrefix(outputFile.getAbsolutePath())+PecanJobData.DECOY_FILE_SUFFIX);
			File featureFile=new File(PecanJobData.getOutputAbsolutePathPrefix(outputFile.getAbsolutePath())+PecanJobData.FEATURE_FILE_SUFFIX);

			try {
				FileLogRecorder logRecorder=new FileLogRecorder(new File(outputFile.getAbsolutePath()+EncyclopediaJobData.LOG_FILE_SUFFIX));
				Logger.addRecorder(logRecorder);
				
				PecanSearchParameters parameters=PecanParameterParser.parseParameters(arguments);
				PecanScoringFactory factory=new PecanOneScoringFactory(parameters, featureFile);
				Logger.logLine("Walnut version "+factory.getVersion());
	
				ArrayList<FastaPeptideEntry> targets;
				if (arguments.containsKey(TARGET_FASTA_TAG)) {
					File targetsFile=new File(arguments.get(TARGET_FASTA_TAG));
					if (SearchParameterParser.getBoolean("-tp", arguments, false)) {
						targets=FastaReader.readPeptideFasta(targetsFile);
					} else {
						targets=new ArrayList<FastaPeptideEntry>();
						ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(targetsFile);
						for (FastaEntryInterface entry : entries) {
							ArrayList<String> peptides=parameters.getEnzyme().digestProtein(entry.getSequence(), parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages(), parameters.getAAConstants().getVariableMods());
							for (String peptide : peptides) {
								FastaPeptideEntry pe=entry.getSubEntry(peptide);
								targets.add(pe);
							}
						}
					}
				} else {
					targets=null;
				}
	
				Logger.logLine("Parameters:");
				Logger.logLine(" "+INPUT_DIA_TAG+" "+diaFile.getAbsolutePath());
				Logger.logLine(" "+BACKGROUND_FASTA_TAG+" "+fastaFile.getAbsolutePath());
				Logger.logLine(" "+TARGET_FASTA_TAG+" "+arguments.get(TARGET_FASTA_TAG));
				Logger.logLine(" "+OUTPUT_RESULT_TAG+" "+outputFile.getAbsolutePath());
				Logger.logLine(parameters.toString());

				runPie(new EmptyProgressIndicator(), Optional.ofNullable(targets), diaFile, fastaFile, featureFile, outputFile, decoyFile, factory);
			} catch (Exception e) {
				Logger.errorLine("Encountered Fatal Error!");
				Logger.errorException(e);
			} finally {
				Logger.close();
			}
		}
	}
	public static void runPie(ProgressIndicator progress, PecanJobData jobData) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		File outputFile=jobData.getOutputFile();
		if (outputFile.exists()&&outputFile.canRead()) {
			try {
				ArrayList<PercolatorPeptide> passingPeptidesFromTSV=PercolatorReader.getPassingPeptidesFromTSV(outputFile, jobData.getParameters().getEffectivePercolatorThreshold(), false);
				ArrayList<ProteinGroup> proteins=ParsimonyProteinGrouper.groupProteins(passingPeptidesFromTSV);
				progress.update("Previously found "+passingPeptidesFromTSV.size()+" peptides ("+proteins.size()+" proteins) identified at "+(jobData.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
				return;
			} catch (Exception e) {
				// problem! so just continue on and overwrite old result
			}
		}
		
		runPie(progress, jobData.getTargetList(), jobData::getDiaFileReader, jobData.getFastaFile(), jobData.getFeatureFile(), jobData.getOutputFile(), jobData.getOutputDecoyFile(), jobData.getTaskFactory());
	}

	static void runPie(ProgressIndicator progress, Optional<ArrayList<FastaPeptideEntry>> targetList, File diaFile, File fastaFile, File featureFile, File outputFile, File decoyFile, PecanScoringFactory taskFactory) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		runPie(
				progress,
				targetList,
				() -> StripeFileGenerator.getFile(diaFile, taskFactory.getParameters()),
				fastaFile,
				featureFile,
				outputFile,
				decoyFile,
				taskFactory
		);
	}

	static void runPie(ProgressIndicator progress, Optional<ArrayList<FastaPeptideEntry>> targetList, Supplier<StripeFileInterface> diaReaderSupplier, File fastaFile, File featureFile, File outputFile, File decoyFile, PecanScoringFactory taskFactory) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		long startTime=System.currentTimeMillis();
		PSMScorer backgroundScorer=taskFactory.getBackgroundScorer();
		PSMPeakScorer pecanScorer=taskFactory.getPecanScorer();
		PecanSearchParameters parameters=taskFactory.getParameters();
		
		int cores=parameters.getNumberOfThreadsUsed();

		Logger.logLine("Converting files...");
		progress.update("Converting files...", Float.MIN_VALUE);

		final StripeFileInterface stripefile = diaReaderSupplier.get();

		Logger.logLine("Processing precursors scans...");
		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));

		PeptideDatabase targets=new PeptideDatabase();
		HashSet<String> backgroundProteome=new HashSet<String>();

		// pecan generates backgrounds using unique fasta peptides, target/decoy sequences, and 2000 random decoys for each window
		// add targets to proteome
		if (targetList.isPresent()) {
			for (FastaPeptideEntry target : targetList.get()) {
				targets.add(target);
				backgroundProteome.add(target.getSequence());
			}
		}

		Logger.logLine("Reading FASTA peptides...");
		// add database to proteome
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(fastaFile);
		for (FastaEntryInterface entry : entries) {
			ArrayList<String> peptides=parameters.getEnzyme().digestProtein(entry.getSequence(), parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages(), parameters.getAAConstants().getVariableMods());
			backgroundProteome.addAll(peptides);

			if (!targetList.isPresent()) {
				// search all peptides in database
				for (String peptide : peptides) {
					FastaPeptideEntry pe=entry.getSubEntry(peptide);
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
			if (!parameters.useTargetWindowCenter()||range.contains(parameters.getTargetWindowCenter())) {
				if (arePeptidesInRange(targets, range, parameters)) {
					ranges.add(range);
				}
			}
		}
		Collections.sort(ranges);
		
		double[] binBoundaries=boundaries.toArray();
		boolean[] useBin=new boolean[binBoundaries.length];
		Arrays.sort(binBoundaries);

		for (Range range : ranges) {
			int index=Arrays.binarySearch(binBoundaries, range.getMiddle());
			if (index>=0) {
				Logger.errorLine("Warning, found window middles that fall on bin boundaries. This implies that the file wasn't demultiplexed correctly!");
			} else {
				index=(-(index+1))-1;
			}
			useBin[index]=true;
		}
		
		Triplet<TDoubleIntHashMap[], ArrayList<String>[], HashSet<String>[]> background=BackgroundGenerator.generateBackground(binBoundaries, useBin, targets, backgroundProteome, parameters);
		TDoubleIntHashMap[] binCounters=background.x;
		ArrayList<String>[] backgroundProteomes=background.y;
		HashSet<String>[] backgroundDecoys=background.z;
		
		BlockingQueue<PeptideScoringResult> resultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
		PeptideScoringResultsConsumer resultsConsumer=taskFactory.getResultsConsumer(resultsQueue, stripefile);
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
			if (index>=0) {
				Logger.errorLine("Warning, found window middles that fall on bin boundaries. This implies that the file wasn't demultiplexed correctly!");
			} else {
				index=(-(index+1))-1;
			}
			TDoubleIntHashMap map=binCounters[index];
			double[] keys=map.keys();
			Arrays.sort(keys);
			ArrayList<String> backgroundProteomeArray=backgroundProteomes[index];
			HashSet<String> backgroundProteomeSet=new HashSet<String>(backgroundProteomeArray);
			
			float dutyCycle=stripefile.getRanges().get(range);
			int scanAveragingMargin=Math.round(parameters.getMinEluteTime()/dutyCycle);
			if (scanAveragingMargin==0) scanAveragingMargin=1;
			
			float maxFragmentationMz=(float)Math.ceil(range.getMiddle()/10.0f)*20.0f+50.0f;
			Range fragmentationRange=new Range(maxFragmentationMz/15f, maxFragmentationMz);
			
			Logger.logLine("Processing "+range+" ("+scanAveragingMargin+")");
			
			List<Stripe> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			Collections.sort(stripes);

			// prepare executor for background
			ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("STRIPE_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
			LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
			ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

			ArrayList<Future<HashMap<LibraryEntry, PeptideScoringResult>>> results=new ArrayList<Future<HashMap<LibraryEntry, PeptideScoringResult>>>();

			int count=0;
			for (String peptide : backgroundDecoys[index]) {
				for (byte charge=parameters.getMinCharge(); charge<=parameters.getMaxCharge(); charge++) {
					double mz=parameters.getAAConstants().getChargedMass(peptide, charge);

					if (range.contains((float)mz)) {
						count++;
						String random=PeptideUtils.getDecoy(peptide, backgroundProteomeSet, parameters);
						AbstractPecanFragmentationModel randmodel=taskFactory.getFragmentationModel(new FastaPeptideEntry(random), parameters.getAAConstants());
						PecanLibraryEntry randentry=randmodel.getPecanSpectrum(charge, keys, map, fragmentationRange, parameters, true);

						ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
						tasks.add(randentry);

						Future<HashMap<LibraryEntry, PeptideScoringResult>> value=executor.submit(new PeptideScoringTask(backgroundScorer, tasks, stripes, precursors, parameters.getAAConstants()));
						results.add(value);
					}
				}
			}
			executor.shutdown();
			while (!executor.isTerminated()) {
				Logger.logLine(workQueue.size()+" background peptides remaining for "+range+"...");
				float finishedFraction=(count-workQueue.size())/(float)count;
				progress.update(baseMessage, baseProgress+baseIncrement*finishedFraction*0.2f);
				Thread.sleep(500);
			}
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

			// prepare executor for peptides
			executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 
			
			
			@SuppressWarnings("unchecked")
			TDoubleObjectHashMap<TDoubleArrayList>[] backgroundScoreMap=new TDoubleObjectHashMap[parameters.getMaxCharge()];
			// background separated out by charge state, TODO assumes positive charge! 
			for (int i=0; i<backgroundScoreMap.length; i++) {
				backgroundScoreMap[i]=new TDoubleObjectHashMap<TDoubleArrayList>();
			}
			
			for (Future<HashMap<LibraryEntry, PeptideScoringResult>> future : results) {
				HashMap<LibraryEntry, PeptideScoringResult> result=future.get();
				for (Entry<LibraryEntry, PeptideScoringResult> resultEntry : result.entrySet()) {
					byte precursorCharge=resultEntry.getKey().getPrecursorCharge();
					
					Pair<double[], double[]> arrays=resultEntry.getValue().getTrace().toArrays();
					double[] x=arrays.x;
					double[] y=arrays.y;
					for (int i=0; i<x.length; i++) {
						TDoubleArrayList list=backgroundScoreMap[precursorCharge-1].get(x[i]);
						if (list==null) {
							list=new TDoubleArrayList();
							backgroundScoreMap[precursorCharge-1].put(x[i], list);
						}
						list.add(y[i]);
					}
				}
			}
			results.clear();
			
			final ArrayList<XYPoint> means=new ArrayList<XYPoint>();

			@SuppressWarnings("unchecked")
			final TDoubleObjectHashMap<XYPoint>[] backgroundScores=new TDoubleObjectHashMap[parameters.getMaxCharge()];
			// background separated out by charge state, TODO assumes positive charge!
			for (int i=0; i<backgroundScores.length; i++) {
				backgroundScores[i]=new TDoubleObjectHashMap<XYPoint>();
			}
			for (int i=0; i<backgroundScores.length; i++) {
				final int chargeIndex=i;
				if (backgroundScoreMap[i].size()>0) {
					backgroundScoreMap[i].forEachEntry(new TDoubleObjectProcedure<TDoubleArrayList>() {
						public boolean execute(double arg0, TDoubleArrayList arg1) {
							double[] values=arg1.toArray();
							double m=General.mean(values);
							double s=General.stdev(values);
							backgroundScores[chargeIndex].put(arg0, new XYPoint(m, s));
							means.add(new XYPoint(arg0, m));
							return true;
						};
					});
				}
			}
			//Charter.launchChart("RT ("+range+" M/Z)", "Fragment Intensity", true, new XYTrace(means, GraphType.line, "Background"));

			count=0;
			for (FastaPeptideEntry peptide : targets) {
				String sequence=peptide.getSequence();
				for (byte charge=parameters.getMinCharge(); charge<=parameters.getMaxCharge(); charge++) {
					double mz=parameters.getAAConstants().getChargedMass(sequence, charge);
					if (range.contains((float)mz)) {
						count++;
						ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
						
						AbstractPecanFragmentationModel model=taskFactory.getFragmentationModel(peptide, parameters.getAAConstants());
						PecanLibraryEntry pecanEntry=model.getPecanSpectrum(charge, keys, map, fragmentationRange, parameters, false);
						tasks.add(pecanEntry);
						
						if (!parameters.isDontRunDecoys()) {
							String smartDecoy=PeptideUtils.getSmartDecoy(sequence, charge, backgroundProteomeSet, parameters);
							FastaPeptideEntry decoyPeptide=new FastaPeptideEntry(peptide.getFilename(), LibraryEntry.DECOY_STRING+peptide.getAccession(), smartDecoy);
							AbstractPecanFragmentationModel revmodel=taskFactory.getFragmentationModel(decoyPeptide, parameters.getAAConstants());
							PecanLibraryEntry reventry=revmodel.getPecanSpectrum(charge, keys, map, fragmentationRange, parameters, true);
							tasks.add(reventry);

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

								String shuffledSequence=PeptideUtils.shuffle(sequence, Float.hashCode(extraDecoys), parameters);
								FastaPeptideEntry shuffledPeptide=new FastaPeptideEntry(peptide.getFilename(), LibraryEntry.SHUFFLE_STRING+peptide.getAccession(), shuffledSequence);
								revmodel=taskFactory.getFragmentationModel(shuffledPeptide, parameters.getAAConstants());
								reventry=revmodel.getPecanSpectrum(charge, keys, map, fragmentationRange, parameters, false);
								tasks.add(reventry);
								
								smartDecoy=PeptideUtils.getSmartDecoy(shuffledSequence, charge, backgroundProteomeSet, parameters);
								decoyPeptide=new FastaPeptideEntry(peptide.getFilename(), LibraryEntry.DECOY_STRING+LibraryEntry.SHUFFLE_STRING+peptide.getAccession(), smartDecoy);
								revmodel=taskFactory.getFragmentationModel(decoyPeptide, parameters.getAAConstants());
								reventry=revmodel.getPecanSpectrum(charge, keys, map, fragmentationRange, parameters, true);
								tasks.add(reventry);
							}
						}

						executor.submit(taskFactory.getScoringTask(pecanScorer, tasks, stripes, backgroundScores, precursors, scanAveragingMargin, resultsQueue));
					}
				}
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
		ArrayList<PercolatorPeptide> passingPeptides=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorVersionNumber(), featureFile, outputFile, decoyFile, parameters.getEffectivePercolatorThreshold());
		ArrayList<ProteinGroup> proteins=ParsimonyProteinGrouper.groupProteins(passingPeptides);
		stripefile.close();
		
		Logger.logLine("Finished analysis! "+resultsConsumer.getNumberProcessed()+" total peaks processed, "+passingPeptides.size()+" peptides ("+proteins.size()+" proteins) identified at 1% FDR ("+(Math.round((System.currentTimeMillis()-startTime)/1000f/6f)/10f)+" minutes)");
		Logger.logLine(""); 
		progress.update(passingPeptides.size()+" peptides identified at "+(parameters.getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
	}

	public static boolean arePeptidesInRange(PeptideDatabase targets, Range range, PecanSearchParameters parameters) {
		// first check to see if we need to process this stripe
		boolean hasPeptides=false;
		outer:for (FastaEntryInterface peptide : targets) {
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
