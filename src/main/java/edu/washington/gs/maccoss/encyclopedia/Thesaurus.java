package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.TreeMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.zip.DataFormatException;

import com.google.common.base.Joiner;
import com.google.common.util.concurrent.ThreadFactoryBuilder;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackground;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.LocalizationDataToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoPermuter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.SaveResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.TeeResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.FileLogRecorder;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.VersioningDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class Thesaurus {

	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		if (arguments.size()==0) {
			SearchGUIMain.runGUI(ProgramType.CASiL);
			
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("Thesaurus Help");
			Logger.timelessLogLine("Thesaurus is software for detecting positional phosphopeptide isomers from DIA data.");
			Logger.timelessLogLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .DIA or .MZML file");
			Logger.timelessLogLine("\t-f\tprotein .FASTA database");
			Logger.timelessLogLine("\t-l\tlibrary .ELIB file");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-o\toutput report file (default: [input file]"+ThesaurusJobData.OUTPUT_FILE_SUFFIX+")");
			
			TreeMap<String, String> defaults=new TreeMap<String, String>(SearchParameterParser.getDefaultParameters());
			int maxWidth=0;
			for (String key : defaults.keySet()) {
				if (key.length()>maxWidth) maxWidth=key.length();
			}
			for (Entry<String, String> entry : defaults.entrySet()) {
				Logger.timelessLogLine("\t"+General.formatCellToWidth(entry.getKey(), maxWidth)+" (default: "+entry.getValue()+")");
			}
			System.exit(1);
			
		} else if (arguments.containsKey("-v")||arguments.containsKey("-version")||arguments.containsKey("--version")) {
			Logger.logLine("Thesaurus version "+ThesaurusOneScoringFactory.version);
			System.exit(1);
			
		} else {
			VersioningDetector.checkVersionCLI(ProgramType.CASiL);
			
			if (!arguments.containsKey(Encyclopedia.INPUT_DIA_TAG)||!arguments.containsKey(Encyclopedia.TARGET_LIBRARY_TAG)||!arguments.containsKey(Encyclopedia.BACKGROUND_FASTA_TAG)) {
				Logger.errorLine("You are required to specify an input file ("+Encyclopedia.INPUT_DIA_TAG+"), a library file ("+Encyclopedia.TARGET_LIBRARY_TAG+"), and a fasta file ("+Encyclopedia.BACKGROUND_FASTA_TAG+")");
				System.exit(1);
			}

			File diaFile=new File(arguments.get(Encyclopedia.INPUT_DIA_TAG));
			File libraryFile=new File(arguments.get(Encyclopedia.TARGET_LIBRARY_TAG));
			File fastaFile=new File(arguments.get(Encyclopedia.BACKGROUND_FASTA_TAG));

			File outputFile;
			if (arguments.containsKey(Encyclopedia.OUTPUT_RESULT_TAG)) {
				outputFile=new File(arguments.get(Encyclopedia.OUTPUT_RESULT_TAG));
			} else {
				outputFile=new File(diaFile.getAbsolutePath()+ThesaurusJobData.OUTPUT_FILE_SUFFIX);
			}

			try {
				FileLogRecorder logRecorder=new FileLogRecorder(new File(outputFile.getAbsolutePath()+ThesaurusJobData.LOG_FILE_SUFFIX));
				Logger.addRecorder(logRecorder);
	
				ThesaurusSearchParameters parameters=ThesaurusSearchParameters.parseParameters(arguments);
				if (!parameters.getLocalizingModification().isPresent()) {
					AminoAcidConstants constants = parameters.getAAConstants();
					String message = getRequiredLocalizationMessage(constants.getLocalizationModifications());
					Logger.errorLine(message);
					System.exit(1);
				}

				Logger.logLine("Setting up localization engine...");
				StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
				PhosphoLocalizer localizer=new PhosphoLocalizer(stripefile, parameters.getLocalizingModification().get(), parameters);
				LibraryScoringFactory factory=new ThesaurusOneScoringFactory(parameters, localizer, new LinkedBlockingQueue<ModificationLocalizationData>());
				
				Logger.logLine("Thesaurus version "+factory.getVersion());
	
				Logger.logLine("Parameters:");
				Logger.logLine(" "+Encyclopedia.INPUT_DIA_TAG+" "+diaFile.getAbsolutePath());
				Logger.logLine(" "+Encyclopedia.TARGET_LIBRARY_TAG+" "+libraryFile.getAbsolutePath());
				Logger.logLine(" "+Encyclopedia.OUTPUT_RESULT_TAG+" "+outputFile.getAbsolutePath());
				Logger.logLine(parameters.toString());

				LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
				ThesaurusJobData job=new ThesaurusJobData(diaFile, library, outputFile, fastaFile, factory);
				runSearch(new EmptyProgressIndicator(), job);
			} catch (Exception e) {
				Logger.errorLine("Encountered Fatal Error!");
				Logger.errorException(e);
			} finally {
				Logger.close();
			}
		}
	}

	public static void runSearch(ProgressIndicator progress, ThesaurusJobData job) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {

		if (job.getPercolatorFiles().hasDataAvailable()) {
			try {
				ArrayList<PercolatorPeptide> passingPeptidesFromTSV=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideOutputFile(), job.getParameters(), false).x;
				
				File elibFile=job.getResultLibrary();
				if (!elibFile.exists()) {
					//job=checkJob(job);
					progress.update("Writing elib result library...");
					Logger.logLine("Writing elib result library...");
					ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
					jobs.add(job);
					SearchToBLIB.convert(progress, jobs, elibFile, false, false);
				}
				Logger.logLine("Previously found "+passingPeptidesFromTSV.size()+" peptides identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR");
				progress.update("Previously found "+passingPeptidesFromTSV.size()+" peptides identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
				//progress.update("Previously found "+passingPeptidesFromTSV.size()+" peptides ("+ParsimonyProteinGrouper.groupProteins(passingPeptidesFromTSV).size()+" proteins) identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);

				return;
			} catch (Exception e) {
				// problem! so just continue on and overwrite old result
				Logger.logLine("Found unexpected exception trying to read old results: ");
				Logger.logException(e);
				Logger.logLine("Just going to go ahead and reprocess this file!");
			}
		}
		
		File diaFile=job.getDiaFile();
		
		Logger.logLine("Converting files...");
		progress.update("Converting files...", Float.MIN_VALUE);
		job=checkJob(job);
		
		SearchParameters parameters=job.getParameters();
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		runSearch(progress, job, stripefile);
		stripefile.close();
	}

	public static ThesaurusJobData checkJob(ThesaurusJobData job) throws IOException, DataFormatException, SQLException {
		if (!(job.getTaskFactory() instanceof ThesaurusOneScoringFactory)) {
			if (!job.getParameters().getLocalizingModification().isPresent()) {
				AminoAcidConstants constants = job.getParameters().getAAConstants();
				String message = getRequiredLocalizationMessage(constants.getLocalizationModifications());
				throw new EncyclopediaException(message);
			}
			ThesaurusSearchParameters searchParameters;
			if (job.getParameters() instanceof ThesaurusSearchParameters) {
				searchParameters=(ThesaurusSearchParameters)job.getParameters();
			} else {
				searchParameters=ThesaurusSearchParameters.convertFromEncyclopeDIA(job.getParameters(), false);
			}

			Logger.logLine("Setting up localization engine...");
			StripeFileInterface stripefile=StripeFileGenerator.getFile(job.getDiaFile(), searchParameters);
			PhosphoLocalizer localizer=new PhosphoLocalizer(stripefile, searchParameters.getLocalizingModification().get(), searchParameters);
			LibraryScoringFactory factory=new ThesaurusOneScoringFactory(searchParameters, localizer, new LinkedBlockingQueue<ModificationLocalizationData>());
			job=job.updateTaskFactory(factory);
		}
		return job;
	}
	
	public static void runSearch(ProgressIndicator progress, ThesaurusJobData job, StripeFileInterface stripefile) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		long startTime=System.currentTimeMillis();
		if (!(job.getTaskFactory() instanceof ThesaurusOneScoringFactory)) {
			throw new EncyclopediaException("Sorry, CASiL requires it's own task factory!");
		}
		
		ThesaurusOneScoringFactory taskFactory=(ThesaurusOneScoringFactory)job.getTaskFactory();
		BlockingQueue<ModificationLocalizationData> localizationQueue=taskFactory.getLocalizationQueue();
		
		SearchParameters parameters=taskFactory.getParameters();

		if (parameters.getNumberOfExtraDecoyLibrariesSearched()>0) {
			Logger.errorLine("Thesaurus does not respect adding extra decoys! Parameter ignored.");
		}
		
		LibraryInterface library=job.getLibrary();
		File featureFile=job.getPercolatorFiles().getInputTSV();
		File localizationFile=job.getLocalizationFile();
		
		int cores=parameters.getNumberOfThreadsUsed();
	
		Logger.logLine("Processing precursors scans...");
		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));
	
		// get targeted ranges
		ArrayList<Range> ranges=new ArrayList<Range>();
		for (Range range : stripefile.getRanges().keySet()) {
			if (!parameters.useTargetWindowCenter()||range.contains(parameters.getTargetWindowCenter())) {
				ranges.add(range);
			}
		}
		Collections.sort(ranges);
	
		PeptideScoringResultsConsumer writeResultsConsumer=taskFactory.getResultsConsumer(featureFile, new LinkedBlockingQueue<PeptideScoringResult>(), stripefile);
		SaveResultsConsumer saveResultsConsumer=new SaveResultsConsumer(new LinkedBlockingQueue<PeptideScoringResult>());
		
		BlockingQueue<PeptideScoringResult> resultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
		TeeResultsConsumer teeConsumer=new TeeResultsConsumer(resultsQueue, writeResultsConsumer, saveResultsConsumer);
		LocalizationDataToTSVConsumer localizationConsumer=new LocalizationDataToTSVConsumer(localizationFile, localizationQueue);
		Thread consumer1Thread=new Thread(teeConsumer);
		Thread consumer2Thread=new Thread(writeResultsConsumer);
		Thread consumer3Thread=new Thread(saveResultsConsumer);
		Thread consumer4Thread=new Thread(localizationConsumer);
		consumer1Thread.start();
		consumer2Thread.start();
		consumer3Thread.start();
		consumer4Thread.start();
	
		// get stripes
		int rangesFinished=0;
		float numberOfTasks=2.0f+ranges.size();
		int totalTasks=0;
		for (Range range : ranges) {
			String baseMessage="Working on "+range+" m/z";
			float baseIncrement=1.0f/numberOfTasks;
			float baseProgress=(1.0f+rangesFinished)/numberOfTasks;
			progress.update(baseMessage, baseProgress);
			
			float dutyCycle=stripefile.getRanges().get(range);
			if (dutyCycle <= 0f) {
				// A stripe with only one scan will get duty cycle
				// of zero. This will only happen in the case of a
				// bad file, or DDA data (where precursor ranges are
				// typically unique). Note that this doesn't guard
				// against (positive) infinity or NaN, but if these
				// values occur it's unclear how to interpret them.
				continue;
			} else if (dutyCycle<=0.1f) {
				dutyCycle=0.1f;
			}
			Logger.logLine("Processing "+range+" m/z, ("+dutyCycle+" second duty cycle)");
			
			ArrayList<FragmentScan> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			Collections.sort(stripes);

			if (stripes.size() < 3) {
				// A stripe with very few scans indicates that either
				// the file is bad, or this is DDA data. Similar to
				// above, we simply skip this stripe.
				continue;
			}

			// prepare executor for background
			ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("STRIPE_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
			LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
			ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 
	
			ArrayList<LibraryEntry> entries=library.getEntries(range, true, parameters.getAAConstants());
			entries=filterEntriesForPTM(entries, parameters.getLocalizingModification(), parameters.getAAConstants());
			LibraryBackgroundInterface background=new LibraryBackground(entries);
			PSMScorer scorer=taskFactory.getLibraryScorer(background);
			
			// keep all entries belonging to the same sequence together
			HashMap<String, ArrayList<LibraryEntry>> entriesBySequence=new HashMap<>();
			for (LibraryEntry entry : entries) {
				String seq=entry.getPeptideSeq();
				ArrayList<LibraryEntry> list=entriesBySequence.get(seq);
				if (list==null) {
					list=new ArrayList<>();
					entriesBySequence.put(seq, list);
				}
				list.add(entry);
			}
			int queueSize=0;
			for (ArrayList<LibraryEntry> entryList : entriesBySequence.values()) {
				queueSize++;
				
				ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
				for (LibraryEntry entry : entryList) {
					tasks.add(entry);
					tasks.add(entry.getDecoy(parameters));
				}
				executor.submit(taskFactory.getScoringTask(scorer, tasks, stripes, range, dutyCycle, precursors, resultsQueue));
				totalTasks++;
			}
			
			executor.shutdown();
			while (!executor.isTerminated()) {
				Logger.logLine(workQueue.size()+" peptides remaining for "+range+"...");
				float finishedFraction=(queueSize-workQueue.size())/(float)queueSize;
				progress.update(baseMessage, baseProgress+baseIncrement*(0.2f+finishedFraction*0.8f));
				Thread.sleep(500);
			}
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
			
			rangesFinished++;
		}
		resultsQueue.put(PeptideScoringResult.POISON_RESULT);
		localizationQueue.put(ModificationLocalizationData.POISON_RESULT);
	
		consumer1Thread.join();
		consumer2Thread.join();
		consumer3Thread.join();
		consumer4Thread.join();
		teeConsumer.close();
		progress.update("Organizing results", (1.0f+rangesFinished)/numberOfTasks);
	
		Logger.logLine("Running Percolator...");
		if (totalTasks==0) {
			throw new EncyclopediaException("No peptides found in library with the specified PTM ("+parameters.getLocalizingModification()+")!");
		}
		Pair<ArrayList<PercolatorPeptide>, RetentionTimeAlignmentInterface> percolatorResults=Encyclopedia.percolatePeptides(progress, job, stripefile, saveResultsConsumer);
		if (parameters.getScoringBreadthType().runRecalibration()) {
			percolatorResults=Encyclopedia.repercolatePeptides(progress, job, stripefile, saveResultsConsumer, percolatorResults.y);
		}
		ArrayList<PercolatorPeptide> passingPeptides=percolatorResults.x;

		localizationConsumer.close();
		
		Logger.logLine("Writing elib result library...");
		try {
			SearchToBLIB.convertElib(progress, job, job.getResultLibrary(), parameters);
		} catch (Exception e) {
			Logger.errorLine("Encountered error creating elib report...");
			Logger.errorException(e);
			throw new EncyclopediaException("Error creating elib report", e);
		}
		
		progress.update("Found "+passingPeptides.size()+" peptides identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
		Logger.logLine("Finished analysis! "+writeResultsConsumer.getNumberProcessed()+" total peptides processed, "+passingPeptides.size()+" peptides identified at "+(parameters.getPercolatorThreshold()*100f)+"% FDR ("+(Math.round((System.currentTimeMillis()-startTime)/1000f/6f)/10f)+" minutes)");
		Logger.logLine(""); 
	}
	
	private static ArrayList<LibraryEntry> filterEntriesForPTM(ArrayList<LibraryEntry> entries, Optional<PeptideModification> optionalptm, AminoAcidConstants aaConstants) {
		if (!optionalptm.isPresent()) return entries;
		
		PeptideModification ptm=optionalptm.get();
		ArrayList<LibraryEntry> pass=new ArrayList<>();
		for (LibraryEntry libraryEntry : entries) {
			if (PhosphoPermuter.hasPTM(libraryEntry.getPeptideModSeq(), ptm, aaConstants)) {
				pass.add(libraryEntry);
			}
		}
		return pass;
	}

	private static String getRequiredLocalizationMessage(Collection<PeptideModification> localizationModifications) {
		String availableLocalizationModifications = Joiner.on(", ").join(
				localizationModifications
						.stream()
						.map(PeptideModification::getShortname)
						.collect(Collectors.toList()));
		return "You are required to specify one localization modification (" + availableLocalizationModifications + ")";
	}

}
