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

import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackground;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.LocalizationDataToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.VariantXCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorrLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorrStripe;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.SimilarPeptideBinner;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.VariantXCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideDatabase;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.parameters.InstrumentSpecificSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.ParsingUtils;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
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

public class VariantXCorDIA {
	public static final String TARGET_FASTA_TAG="-t";
	public static final String TARGET_LIBRARY_TAG="-l";
	public static final String OUTPUT_RESULT_TAG="-o";
	public static final String INPUT_DIA_TAG="-i";
	public static final String BACKGROUND_FASTA_TAG="-f";

	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		arguments=InstrumentSpecificSearchParameters.checkParameters(arguments);
		
		if (arguments.size()==0) {
			SearchGUIMain.runGUI(ProgramType.XCorDIA);
			
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("XCorDIA Help");
			Logger.timelessLogLine("XCorDIA is a FASTA database search engine for DIA data.");
			Logger.timelessLogLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .DIA or .MZML file");
			Logger.timelessLogLine("\t-f\tbackground FASTA file");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-t\ttarget FASTA file (default: background FASTA file, not used if library specified)");
			Logger.timelessLogLine("\t-l\ttarget Library file (default: search target FASTA)");
			Logger.timelessLogLine("\t-tp\ttrue/false target FASTA file contains peptides (default: false)"); 
			Logger.timelessLogLine("\t-o\toutput report file (default: [input file]"+XCorDIAJobData.OUTPUT_FILE_SUFFIX+")");

			TreeMap<String, String> defaults=new TreeMap<String, String>(PecanParameterParser.getDefaultParameters());
			int maxWidth=0;
			for (String key : defaults.keySet()) {
				if (key.length()>maxWidth) maxWidth=key.length();
			}
			for (Entry<String, String> entry : defaults.entrySet()) {
				Logger.timelessLogLine("\t"+General.formatCellToWidth(entry.getKey(), maxWidth)+" (default: "+entry.getValue()+")");
			}

			Logger.timelessLogLine("\t"+Encyclopedia.QUIET_MODE_ARG+"\tsuppress log output to stdout/stderr");

			System.exit(1);
			
		} else if (arguments.containsKey("-v")||arguments.containsKey("-version")||arguments.containsKey("--version")) {
			Logger.logLine("XCorDIA version "+ProgramType.getGlobalVersion().toString());
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
				outputFile=new File(diaFile.getAbsolutePath()+XCorDIAJobData.OUTPUT_FILE_SUFFIX);
			}

			try {
				if (arguments.containsKey(Encyclopedia.QUIET_MODE_ARG)) {
					Logger.PRINT_TO_SCREEN = false;
				}
				FileLogRecorder logRecorder=new FileLogRecorder(new File(outputFile.getAbsolutePath()+EncyclopediaJobData.LOG_FILE_SUFFIX));
				Logger.addRecorder(logRecorder);
				
				PecanSearchParameters parameters=PecanParameterParser.parseParameters(arguments);
				XCorDIAOneScoringFactory factory=new XCorDIAOneScoringFactory(parameters);
				Logger.logLine("XCorDIA version "+ProgramType.getGlobalVersion().toString());
	
				ArrayList<FastaPeptideEntry> targets;
				if (arguments.containsKey(TARGET_FASTA_TAG)) {
					File targetsFile=new File(arguments.get(TARGET_FASTA_TAG));
					if (ParsingUtils.getBoolean("-tp", arguments, false)) {
						targets=FastaReader.readPeptideFasta(targetsFile, parameters);
					} else {
						targets=new ArrayList<FastaPeptideEntry>();
						ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(targetsFile, parameters);
						for (FastaEntryInterface entry : entries) {
							ArrayList<FastaPeptideEntry> peptides=parameters.getEnzyme().digestProtein(entry, parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages(), parameters.getAAConstants(), parameters.isRequireVariableMods());
							for (FastaPeptideEntry peptide : peptides) {
								targets.add(peptide);
							}
						}
					}
				} else {
					targets=null;
				}

				LibraryInterface library;
				if (arguments.containsKey(TARGET_LIBRARY_TAG)) {
					library=BlibToLibraryConverter.getFile(new File(arguments.get(TARGET_LIBRARY_TAG)), fastaFile, parameters);
				} else {
					library=null;
				}
				
				VariantXCorDIAJobData jobData=new VariantXCorDIAJobData(Optional.ofNullable(targets), Optional.ofNullable(library), diaFile, fastaFile, outputFile, factory);

				runPie(new EmptyProgressIndicator(), jobData);
			} catch (Exception e) {
				Logger.errorLine("Encountered Fatal Error!");
				Logger.errorException(e);
			} finally {
				Logger.close();
			}
		}
	}

	public static VariantXCorDIAJobData checkJob(VariantXCorDIAJobData job) throws IOException, DataFormatException, SQLException {
		if (!(job.getTaskFactory() instanceof VariantXCorDIAOneScoringFactory)) {
			Logger.logLine("Setting up localization engine...");
			PecanSearchParameters parameters=job.getTaskFactory().getPecanParameters();
			StripeFileInterface stripefile=job.getDiaFileReader();
			BackgroundFrequencyInterface background=BackgroundFrequencyCalculator.generateBackground(stripefile);
			//BackgroundFrequencyInterface background=new UnitBackgroundFrequencyCalculator(0.01f); //FIXME REMOVE, TESTING ONLY
			VariantXCorDIAOneScoringFactory factory=new VariantXCorDIAOneScoringFactory(parameters, background, new LinkedBlockingQueue<ModificationLocalizationData>());
			job=job.updateTaskFactory(factory);
		}
		return job;
	}
	
	public static void runPie(ProgressIndicator progress, VariantXCorDIAJobData jobData) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		if (jobData.getPercolatorFiles().hasDataAvailable()) {
			try {
				ArrayList<PercolatorPeptide> passingPeptidesFromTSV=PercolatorReader.getPassingPeptidesFromTSV(jobData.getPercolatorFiles().getPeptideOutputFile(), jobData.getParameters().getEffectivePercolatorThreshold(), jobData.getParameters().getAAConstants(), false).x;
				File elibFile=jobData.getResultLibrary();
				if (!elibFile.exists()) {
					//job=checkJob(job);
					progress.update("Writing elib result library...");
					Logger.logLine("Writing elib result library...");
					ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
					jobs.add(jobData);
					SearchToBLIB.convert(progress, jobs, elibFile, false, false);
				}
				progress.update("Previously found "+passingPeptidesFromTSV.size()+" peptides identified at "+(jobData.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
				return;
			} catch (Exception e) {
				Logger.errorException(e);
				Logger.errorLine("Found error reading old results, discarding and reprocessing!");
				// problem! so just continue on and overwrite old result
			}
		}
		jobData=checkJob(jobData);
		
		long startTime=System.currentTimeMillis();
		final PecanSearchParameters parameters=jobData.getTaskFactory().getPecanParameters();
		
		int cores=parameters.getNumberOfThreadsUsed();
		
		Logger.logLine("Using "+jobData.getTaskFactory().getName());
		Logger.logLine("Input File: "+jobData.getOriginalDiaFileName());
		if (jobData.getLibrary().isPresent()) {
			Logger.logLine("Library File: "+jobData.getLibrary().get().getName());
		} else {
			Logger.logLine("Library File Missing.");
		}
		Logger.logLine("FASTA File: "+jobData.getFastaFile().getName());
		Logger.logLine("Result File: "+jobData.getResultLibrary().getName());
		Logger.logLine("Parameters:");
		Logger.logLine(jobData.getParameters().toString());

		Logger.logLine("Converting files...");
		progress.update("Converting files...", Float.MIN_VALUE);

		final StripeFileInterface stripefile = jobData.getDiaFileReader();

		Logger.logLine("Processing precursors scans...");
		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));

		PeptideDatabase targets=new PeptideDatabase();
		HashSet<String> backgroundProteome=new HashSet<String>();

		// xcordia generates backgrounds using unique fasta peptides
		if (jobData.getTargetList().isPresent()) {
			for (FastaPeptideEntry target : jobData.getTargetList().get()) {
				targets.add(target);
				backgroundProteome.add(target.getSequence());
			}
		}
		LibraryBackgroundInterface background=new LibraryBackground(backgroundProteome, parameters);
		XCorDIAOneScorer xcordiaScorer=(XCorDIAOneScorer)jobData.getTaskFactory().getLibraryScorer(background);
		
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

		Logger.logLine("Reading FASTA peptides...");
		// add database to proteome
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(jobData.getFastaFile(), parameters);
		Logger.logLine("Found "+entries.size()+" total proteins...");
		
		for (FastaEntryInterface entry : entries) {
			ArrayList<FastaPeptideEntry> peptides=parameters.getEnzyme().digestProtein(entry, parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages(), parameters.getAAConstants(), parameters.isRequireVariableMods());

			for (FastaPeptideEntry peptide : peptides) {
				backgroundProteome.add(peptide.getSequence());
			}

			if (!jobData.getTargetList().isPresent()) {
				for (FastaPeptideEntry peptide : peptides) {
					targets.add(peptide);
				}
			}
		}
		Logger.logLine("Found "+targets.size()+" total peptides from "+entries.size()+" proteins...");
		
		// get targeted ranges
		float minMz=Float.MAX_VALUE;
		float maxMz=-Float.MAX_VALUE;
		TDoubleHashSet boundaries=new TDoubleHashSet();
		ArrayList<Range> ranges=new ArrayList<Range>();
		for (Range range : stripefile.getRanges().keySet()) {
			if (range.getStart()<minMz) minMz=range.getStart();
			if (range.getStop()>maxMz) maxMz=range.getStop();
			boundaries.add(range.getStart());
			boundaries.add(range.getStop());
			if (!parameters.useTargetWindowCenter()||range.contains(parameters.getTargetWindowCenter())) {
				if (XCorDIA.areAnyPeptidesInRange(targets, range, parameters)) {
					ranges.add(range);
				}
			}
		}

		for (FastaPeptideEntry peptide : targets.getPeptides()) {
			String sequence=peptide.getSequence();

			byte expectedCharge=PeptideUtils.getExpectedChargeState(sequence);
			byte minCharge=(byte)Math.max(parameters.getMinCharge(), expectedCharge-1);
			byte maxCharge=(byte)Math.min(parameters.getMaxCharge(), expectedCharge+1);
				
			boolean keep=false;
			for (byte charge=minCharge; charge<=maxCharge; charge++) {
				double mz=parameters.getAAConstants().getChargedMass(sequence, charge);
				if (mz>minMz&&mz<maxMz) {
					keep=true;
					break;
				}
			}
			if (!keep) {
				targets.remove(sequence);
			}
		}
		Collections.sort(ranges);
		Logger.logLine("Found "+ranges.size()+" matching ranges, trimmed targets down to "+targets.size()+" peptides...");
		
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
		
		// localization queues
		BlockingQueue<ModificationLocalizationData> localizationQueue=((VariantXCorDIAOneScoringFactory)jobData.getTaskFactory()).getLocalizationQueue();
		File localizationFile=jobData.getLocalizationFile();
		
		BlockingQueue<AbstractScoringResult> resultsQueue=new LinkedBlockingQueue<AbstractScoringResult>();
		PeptideScoringResultsConsumer resultsConsumer=jobData.getTaskFactory().getResultsConsumer(jobData.getPercolatorFiles().getInputTSV(), resultsQueue, stripefile, jobData.getLibrary().orElse(null));
		LocalizationDataToTSVConsumer localizationConsumer=new LocalizationDataToTSVConsumer(localizationFile, localizationQueue);
		Thread consumerThread=new Thread(resultsConsumer);
		consumerThread.start();
		Thread consumer4Thread=new Thread(localizationConsumer);
		consumer4Thread.start();
		
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
			
			float dutyCycle=stripefile.getRanges().get(range).getAverageDutyCycle();
			if (dutyCycle <= 0f) {
				// A stripe with only one scan will get duty cycle
				// of zero. This will only happen in the case of a
				// bad file, or DDA data (where precursor ranges are
				// typically unique). Note that this doesn't guard
				// against (positive) infinity or NaN, but if these
				// values occur it's unclear how to interpret them.
				continue;
			}

			Logger.logLine("Processing "+range+" ("+dutyCycle+")");

			// prepare executor
			ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("STRIPE_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
			LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
			ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

			// set up xcorr
			ArrayList<FragmentScan> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);

			if (stripes.size() < 3) {
				// A stripe with very few scans indicates that either
				// the file is bad, or this is DDA data. Similar to
				// above, we simply skip this stripe.
				continue;
			}

			Logger.logLine("Starting XCorr background calculations for "+stripes.size()+" spectra between "+range+"...");
			final Vector<FragmentScan> tempStripes=new Vector<FragmentScan>();
			for (final FragmentScan stripe : stripes) {
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
			
			HashSet<LibraryEntry> allPeptidesInWindow=new HashSet<>();
			if (jobData.getLibrary().isPresent()) {
				HashSet<FastaPeptideEntry> allPeptideSequencesInWindow=XCorDIA.getPeptidesInRange(parameters, targets, range);
				HashSet<String> sequences=new HashSet<>();
				for (FastaPeptideEntry entry : allPeptideSequencesInWindow) {
					sequences.add(entry.getSequenceWithModsStripped());
				}
				
				ArrayList<LibraryEntry> libEntries=jobData.getLibrary().get().getEntries(range, true, parameters.getAAConstants());
				Logger.logLine("Extracting spectra from library, considered "+libEntries.size()+" spectra...");
				
				for (LibraryEntry libraryEntry : libEntries) {
					if (sequences.contains(libraryEntry.getPeptideSeq())) {
						XCorrLibraryEntry entry=XCorrLibraryEntry.generateEntry(false, libraryEntry, libraryEntry.getPrecursorCharge(), parameters);
						allPeptidesInWindow.add(entry);
					}
				}
			} else {
				HashSet<FastaPeptideEntry> allPeptideSequencesInWindow=XCorDIA.getPeptidesInRange(parameters, targets, range);
				Logger.logLine("Extracting spectra from PEFF, considered "+allPeptideSequencesInWindow.size()+" peptides...");
				
				for (FastaPeptideEntry peptide : allPeptideSequencesInWindow) {
					String sequence=peptide.getSequence();
	
					byte expectedCharge=PeptideUtils.getExpectedChargeState(sequence);
					byte minCharge=(byte)Math.max(parameters.getMinCharge(), expectedCharge-1);
					byte maxCharge=(byte)Math.min(parameters.getMaxCharge(), expectedCharge+1);
					
					for (byte charge=minCharge; charge<=maxCharge; charge++) {
						double mz=parameters.getAAConstants().getChargedMass(sequence, charge);
						if (range.contains((float)mz)) {
							XCorrLibraryEntry entry=XCorrLibraryEntry.generateEntry(false, peptide, charge, parameters);
							allPeptidesInWindow.add(entry);
						}
					}
				}
			}
			
			

			SimilarPeptideBinner binner=new SimilarPeptideBinner();
			ArrayList<ArrayList<LibraryEntry>> bins=binner.binEntries(allPeptidesInWindow);
			Logger.logLine("Found "+allPeptidesInWindow.size()+" peptides in "+bins.size()+" peptide groups...");

			int count=0;

			for (ArrayList<LibraryEntry> bin : bins) {
				ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
				count++;
				for (LibraryEntry entry : bin) {
					tasks.add(entry);
					
					if (!parameters.isDontRunDecoys()) {
						LibraryEntry decoy=((XCorrLibraryEntry)entry).getDecoy(parameters);
						tasks.add(decoy);

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

							LibraryEntry shuffle=((XCorrLibraryEntry)entry).getShuffle(parameters, Float.hashCode(extraDecoys), false);
							tasks.add(shuffle);
							tasks.add(((XCorrLibraryEntry)shuffle).getDecoy(parameters));
						}
					}
				}
				executor.submit(jobData.getTaskFactory().getScoringTask(xcordiaScorer, tasks, stripes, range, dutyCycle, precursors, resultsQueue));
			}
			executor.shutdown();
			while (!executor.isTerminated()) {
				Logger.logLine(workQueue.size()+" peptide groups remaining for "+range+"...");
				float finishedFraction=(count-workQueue.size())/(float)count;
				progress.update(baseMessage, baseProgress+baseIncrement*(0.2f+finishedFraction*0.8f));
				Thread.sleep(500);
			}
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
			
			rangesFinished++;
		}
		resultsQueue.put(AbstractScoringResult.POISON_RESULT);
		localizationQueue.put(ModificationLocalizationData.POISON_RESULT);

		consumerThread.join();
		consumer4Thread.join();
		resultsConsumer.close();
		localizationConsumer.close();
		
		Logger.logLine("Finished generating feature file, analyzed "+resultsConsumer.getNumberProcessed()+" peptides.");

		progress.update("Running Percolator", (1.0f+rangesFinished)/numberOfTasks);
		ArrayList<PercolatorPeptide> passingPeptides=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorVersionNumber(), jobData.getPercolatorFiles(), parameters.getEffectivePercolatorThreshold(), parameters.getAAConstants(), 1).x;
		
		Logger.logLine("Writing elib result library...");
		SearchToBLIB.convertElib(progress, jobData, jobData.getResultLibrary(), parameters, parameters.isIntegratePrecursors());
		stripefile.close();
		
		Logger.logLine("Finished analysis! "+resultsConsumer.getNumberProcessed()+" total peaks processed, "+passingPeptides.size()+" peptides identified at 1% FDR ("+(Math.round((System.currentTimeMillis()-startTime)/1000f/6f)/10f)+" minutes)");
		Logger.logLine(""); 
		progress.update(passingPeptides.size()+" peptides identified at "+(parameters.getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
	}
}
