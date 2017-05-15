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
import java.util.Optional;
import java.util.TreeMap;
import java.util.Vector;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ParsimonyProteinGrouper;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackground;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorrLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorrStripe;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideDatabase;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.FileLogRecorder;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ThreadableTask;
import gnu.trove.set.hash.TDoubleHashSet;

public class XCorDIA {
	public static final String TARGET_FASTA_TAG="-t";
	public static final String OUTPUT_RESULT_TAG="-o";
	public static final String INPUT_DIA_TAG="-i";
	public static final String BACKGROUND_FASTA_TAG="-f";

	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		if (arguments.size()==0) {
			SearchGUIMain.runGUI(true);
			
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("XCorDIA Help");
			Logger.timelessLogLine("XCorDIA is a FASTA database search engine for DIA data.");
			Logger.timelessLogLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .DIA or .MZML file");
			Logger.timelessLogLine("\t-f\tbackground FASTA file");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-t\ttarget FASTA file (default: background FASTA file)");
			Logger.timelessLogLine("\t-tp\ttrue/false target FASTA file contains peptides (default: false)"); 
			Logger.timelessLogLine("\t-o\toutput report file (default: [input file].xcordia.txt)");
			
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
			Logger.logLine("XCorDIA version "+XCorDIAOneScoringFactory.version);
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
				outputFile=new File(diaFile.getAbsolutePath()+".xcordia.txt");
			}

			File featureFile=new File(outputFile.getAbsolutePath()+".features.txt");

			try {
				FileLogRecorder logRecorder=new FileLogRecorder(new File(outputFile.getAbsolutePath()+EncyclopediaJobData.LOG_FILE_SUFFIX));
				Logger.addRecorder(logRecorder);
				
				PecanSearchParameters parameters=PecanParameterParser.parseParameters(arguments);
				XCorDIAOneScoringFactory factory=new XCorDIAOneScoringFactory(parameters);
				Logger.logLine("XCorDIA version "+factory.getVersion());
	
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

				runPie(new EmptyProgressIndicator(), Optional.ofNullable(targets), diaFile, fastaFile, featureFile, outputFile, factory);
			} catch (Exception e) {
				Logger.errorLine("Encountered Fatal Error!");
				Logger.errorException(e);
			} finally {
				Logger.close();
			}
		}
	}
	public static void runPie(ProgressIndicator progress, XCorDIAJobData jobData) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		File outputFile=jobData.getOutputFile();
		if (outputFile.exists()&&outputFile.canRead()) {
			try {
				ArrayList<PercolatorPeptide> passingPeptidesFromTSV=PercolatorReader.getPassingPeptidesFromTSV(outputFile, jobData.getParameters().getEffectivePercolatorThreshold());
				ArrayList<ProteinGroup> proteins=ParsimonyProteinGrouper.groupProteins(passingPeptidesFromTSV);
				progress.update("Previously found "+passingPeptidesFromTSV.size()+" peptides ("+proteins.size()+" proteins) identified at "+(jobData.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
				return;
			} catch (Exception e) {
				// problem! so just continue on and overwrite old result
			}
		}
		
		runPie(progress, jobData.getTargetList(), jobData.getDiaFile(), jobData.getFastaFile(), jobData.getFeatureFile(), jobData.getOutputFile(), jobData.getTaskFactory());
	}
		
	static void runPie(ProgressIndicator progress, Optional<ArrayList<FastaPeptideEntry>> targetList, File diaFile, File fastaFile, File featureFile, File outputFile, XCorDIAOneScoringFactory taskFactory) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		long startTime=System.currentTimeMillis();
		final PecanSearchParameters parameters=taskFactory.getPecanParameters();
		
		int cores=parameters.getNumberOfThreadsUsed();

		Logger.logLine("Converting files...");
		progress.update("Converting files...", Float.MIN_VALUE);
		
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);

		Logger.logLine("Processing precursors scans...");
		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));

		PeptideDatabase targets=new PeptideDatabase();
		HashSet<String> backgroundProteome=new HashSet<String>();

		// xcordia generates backgrounds using unique fasta peptides, target/decoy sequences, and 2000 random decoys for each window
		// add targets to proteome
		if (targetList.isPresent()) {
			for (FastaPeptideEntry target : targetList.get()) {
				targets.add(target);
				backgroundProteome.add(target.getSequence());
			}
		}
		LibraryBackgroundInterface background=new LibraryBackground(backgroundProteome, parameters);
		XCorDIAOneScorer xcordiaScorer=(XCorDIAOneScorer)taskFactory.getLibraryScorer(background);

		float minMzRange=Float.MAX_VALUE;
		float maxMzRange=-Float.MAX_VALUE;
		for (Range range : stripefile.getRanges().keySet()) {
			if (minMzRange>range.getStart()) {
				minMzRange=range.getStart();
			}
			if (maxMzRange<range.getStop()) {
				maxMzRange=range.getStop();
			}
		}
//		Range entireRange=new Range(minMzRange, maxMzRange);

		Logger.logLine("Reading FASTA peptides...");
		// add database to proteome
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(fastaFile);
		for (FastaEntryInterface entry : entries) {
			ArrayList<String> peptides=parameters.getEnzyme().digestProtein(entry.getSequence(), parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages(), parameters.getAAConstants().getVariableMods());
			backgroundProteome.addAll(peptides);

			if (!targetList.isPresent()) {
//				if (false) { // FIXME
//					ArrayList<ScoredObject<String>> scoredPeptides=new ArrayList<ScoredObject<String>>();
//					for (String peptide : peptides) {
//						byte expectedCharge=PeptideUtils.getExpectedChargeState(peptide);
//						double mz=parameters.getAAConstants().getChargedMass(peptide, expectedCharge);
//						if (entireRange.contains(mz)) {
//							double score=Prego.getNetwork().getScore(peptide);
//							scoredPeptides.add(new ScoredObject<String>((float)score, peptide));
//						}
//					}
//
//					Collections.sort(scoredPeptides);
//					int count=0;
//					for (int i=scoredPeptides.size()-1; i>=0; i--) {
//						count++;
//						FastaPeptideEntry pe=entry.getSubEntry(scoredPeptides.get(i).y);
//						targets.add(pe);
//						if (count>=10) break;
//					}
//				} else {
					// search all peptides in database
					for (String peptide : peptides) {
						FastaPeptideEntry pe=entry.getSubEntry(peptide);
						targets.add(pe);
					}
//				}
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
		
		BlockingQueue<PeptideScoringResult> resultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
		PeptideScoringResultsConsumer resultsConsumer=taskFactory.getResultsConsumer(featureFile, resultsQueue, diaFile);
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
			
			float dutyCycle=stripefile.getRanges().get(range);
			Logger.logLine("Processing "+range+" ("+dutyCycle+")");

			// prepare executor
			ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("STRIPE_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
			LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
			ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

			// set up xcorr
			ArrayList<Stripe> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);

			Logger.logLine("Starting XCorr background calculations for "+stripes.size()+" spectra between "+range+"...");
			final Vector<Stripe> tempStripes=new Vector<Stripe>();
			for (final Stripe stripe : stripes) {
				executor.submit(new ThreadableTask<Nothing>() {
					@Override
					public String getTaskName() {
						return "XCorr background calculation";
					}
					@Override
					protected Nothing process() {
						tempStripes.add(new XCorrStripe(stripe, parameters));
						return Nothing.NOTHING;
					}
				});
			}
			executor.shutdown();
			while (!executor.isTerminated()) {
				Logger.logLine("Processing XCorr background calculations for "+workQueue.size()+" remaining spectra between "+range+"...");
				float finishedFraction=(stripes.size()-workQueue.size())/(float)stripes.size();
				progress.update(baseMessage, baseProgress+baseIncrement*(finishedFraction*0.2f));
				Thread.sleep(500);
			}
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
			
			stripes.clear();
			stripes.addAll(tempStripes);
			Collections.sort(stripes);
			tempStripes.clear();
			Logger.logLine("Finished XCorr background analysis for "+range+", starting queueing peptides...");
			
			executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 
			
			int count=0;
			for (FastaPeptideEntry peptide : targets) {
				String sequence=peptide.getSequence();

				byte minCharge;
				byte maxCharge;
				byte expectedCharge=PeptideUtils.getExpectedChargeState(sequence);
				//if (false&&!targetList.isPresent()) { //FIXME
				//	minCharge=expectedCharge;
				//	maxCharge=expectedCharge;
				//} else {
					minCharge=(byte)Math.max(parameters.getMinCharge(), expectedCharge-1);
					maxCharge=(byte)Math.min(parameters.getMaxCharge(), expectedCharge+1);
				//}
				for (byte charge=minCharge; charge<=maxCharge; charge++) {
					double mz=parameters.getAAConstants().getChargedMass(sequence, charge);
					if (range.contains((float)mz)) {
						count++;
						ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
						
						XCorrLibraryEntry entry=XCorrLibraryEntry.generateEntry(false, peptide.getFilename(), peptide.getAccessions(), charge, peptide.getSequence(), parameters);
						tasks.add(entry);
						
						if (!parameters.isDontRunDecoys()) {
							String smartDecoy=PeptideUtils.reverse(sequence, parameters);
							FastaPeptideEntry decoyPeptide=new FastaPeptideEntry(peptide.getFilename(), LibraryEntry.DECOY_STRING+peptide.getAccession(), smartDecoy);
							XCorrLibraryEntry reventry=XCorrLibraryEntry.generateEntry(true, decoyPeptide.getFilename(), decoyPeptide.getAccessions(), charge, decoyPeptide.getSequence(), parameters);
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
								reventry=XCorrLibraryEntry.generateEntry(true, shuffledPeptide.getFilename(), shuffledPeptide.getAccessions(), charge, shuffledPeptide.getSequence(), parameters);
								tasks.add(reventry);
								
								smartDecoy=PeptideUtils.reverse(shuffledSequence, parameters);
								decoyPeptide=new FastaPeptideEntry(peptide.getFilename(), LibraryEntry.DECOY_STRING+LibraryEntry.SHUFFLE_STRING+peptide.getAccession(), smartDecoy);
								reventry=XCorrLibraryEntry.generateEntry(true, decoyPeptide.getFilename(), decoyPeptide.getAccessions(), charge, decoyPeptide.getSequence(), parameters);
								tasks.add(reventry);
							}
						}

						executor.submit(taskFactory.getScoringTask(xcordiaScorer, tasks, stripes, dutyCycle, precursors, resultsQueue));
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
		Logger.logLine("Finished generating feature file, analyzed "+resultsConsumer.getNumberProcessed()+" peptides.");

		progress.update("Running Percolator", (1.0f+rangesFinished)/numberOfTasks);
		ArrayList<PercolatorPeptide> passingPeptides=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorVersionNumber(), featureFile, outputFile, parameters.getEffectivePercolatorThreshold());
		stripefile.close();
		
		/*if (false&&!targetList.isPresent()) { //FIXME
			HashSet<String> accessions=new HashSet<String>();
			for (PercolatorPeptide peptide : passingPeptides) {
				accessions.addAll(PSMData.stringToAccessions(peptide.getProteinIDs()));
			}
			
			ArrayList<FastaPeptideEntry> newTargets=new ArrayList<FastaPeptideEntry>();
			for (FastaEntryInterface entry : entries) {
				if (accessions.contains(entry.getAccession())) {
					ArrayList<String> peptides=parameters.getEnzyme().digestProtein(entry.getSequence(), parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages());
					for (String peptide : peptides) {
						FastaPeptideEntry pe=entry.getSubEntry(peptide);
						newTargets.add(pe);
					}
				}
			}
			
			runPie(progress, Optional.of(newTargets), diaFile, fastaFile, featureFile, outputFile, taskFactory);
			
		} else {*/
			ArrayList<ProteinGroup> proteins=ParsimonyProteinGrouper.groupProteins(passingPeptides);
			Logger.logLine("Finished analysis! "+resultsConsumer.getNumberProcessed()+" total peaks processed, "+passingPeptides.size()+" peptides ("+proteins.size()+" proteins) identified at 1% FDR ("+(Math.round((System.currentTimeMillis()-startTime)/1000f/6f)/10f)+" minutes)");
			Logger.logLine(""); 
			progress.update(passingPeptides.size()+" peptides identified at "+(parameters.getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
		//}
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
