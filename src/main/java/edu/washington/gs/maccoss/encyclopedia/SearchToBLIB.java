package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.Lists;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ParsimonyProteinGrouper;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.AlternatePeakLocationInferrer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrerInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.LocalizationDataToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.precursor.DDAPrecursorIntegrator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.LibraryReportExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.PeptideQuantExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.VariantXCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
import edu.washington.gs.maccoss.encyclopedia.filereaders.*;
import edu.washington.gs.maccoss.encyclopedia.utils.*;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableConcatenator;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.*;
import java.util.Map.Entry;
import java.util.zip.DataFormatException;

public class SearchToBLIB {
	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		if (arguments.size()==0) {
			SearchGUIMain.runGUI(ProgramType.EncyclopeDIA);
			
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("SearchToLIB Help");
			Logger.timelessLogLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.timelessLogLine("Other Programs: ");
			Logger.timelessLogLine("\t-pecan\trun Pecanpie export (use -pecan -h for Pecan help)");
			Logger.timelessLogLine("\t-xcordia\trun XCorDIA export (use -xcordia -h for XCorDIA help)");
			Logger.timelessLogLine("\t-phospho\trun phospho localization export (use -phospho -h for localization help)");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .DIA or .MZML file or directory");
			Logger.timelessLogLine("\t-o\toutput library .ELIB file");
			Logger.timelessLogLine("\t-a\talign between files (default=true)");
			Logger.timelessLogLine("\t-blib\twrite .BLIB instead of .ELIB (default=false)");
			Logger.timelessLogLine("Potentially Required Parameters: ");
			Logger.timelessLogLine("\t-l\toriginal searched library .DLIB or .ELIB file (required by EncyclopeDIA Export)");
			Logger.timelessLogLine("\t-f\toriginal fasta file (required by Pecan/XCorDIA Export)");
			Logger.timelessLogLine("\t-t\toriginal target file (optional for Pecan/XCorDIA Export)");

			Logger.timelessLogLine("Other Parameters: ");
			TreeMap<String, String> defaults=new TreeMap<String, String>(SearchParameterParser.getExportParameters());
			int maxWidth=0;
			for (String key : defaults.keySet()) {
				if (key.length()>maxWidth) maxWidth=key.length();
			}
			for (Entry<String, String> entry : defaults.entrySet()) {
				Logger.timelessLogLine("\t"+General.formatCellToWidth(entry.getKey(), maxWidth)+" (default: "+entry.getValue()+")");
			}
			System.exit(1);
			
		} else if (arguments.containsKey("-v")||arguments.containsKey("-version")||arguments.containsKey("--version")) {
			Logger.logLine("EncyclopeDIA SearchToLIB version "+ProgramType.getGlobalVersion().toString());
			System.exit(1);
			
		} else {
			if (arguments.containsKey("-pecan")||arguments.containsKey("-walnut")) {
				VersioningDetector.checkVersionCLI(ProgramType.PecanPie);
				convertPecan(arguments);
			} else if (arguments.containsKey("-xcordia")) {
				VersioningDetector.checkVersionCLI(ProgramType.XCorDIA);
				convertXCorDIA(arguments);
			} else {
				VersioningDetector.checkVersionCLI(ProgramType.EncyclopeDIA);
				convertEncyclopedia(arguments);
			}
		}
	}

	public static void convertXCorDIA(HashMap<String, String> arguments) {
		if (!arguments.containsKey("-i")||!arguments.containsKey("-f")||!arguments.containsKey("-o")) {
			Logger.errorLine("You are required to specify an input file or directory (-i), an input fasta file (-f) and an output library file (-o)");
			System.exit(1);
		}

		File diaFile=new File(arguments.get("-i"));
		File fastaFile=new File(arguments.get("-f"));
		File outputFile=new File(arguments.get("-o"));
		boolean alignBetweenFiles=ParsingUtils.getBoolean("-a", arguments, true);
		boolean writeBlib=ParsingUtils.getBoolean("-blib", arguments, false);

		final OutputFormat outputFormat = writeBlib ? OutputFormat.BLIB : OutputFormat.ELIB;

		PecanSearchParameters parameters=PecanParameterParser.parseParameters(arguments);
		XCorDIAOneScoringFactory factory=new XCorDIAOneScoringFactory(parameters);
		Logger.timelessLogLine("SearchToLIB XCorDIA version "+ProgramType.getGlobalVersion().toString());

		Logger.timelessLogLine("Parameters:");
		Logger.timelessLogLine(" -i "+diaFile.getAbsolutePath());
		Logger.timelessLogLine(" -f "+fastaFile.getAbsolutePath());
		Logger.timelessLogLine(" -o "+outputFile.getAbsolutePath());
		Logger.timelessLogLine(" -a "+alignBetweenFiles);
		Logger.timelessLogLine(" -blib "+writeBlib);
		Logger.timelessLogLine(parameters.toString());

		try {
			ArrayList<FastaPeptideEntry> targets;
			if (arguments.containsKey(XCorDIA.TARGET_FASTA_TAG)) {
				targets=FastaReader.readPeptideFasta(new File(arguments.get(XCorDIA.TARGET_FASTA_TAG)), parameters);
			} else {
				targets=null;
			}
			LibraryInterface library;
			if (arguments.containsKey("-l")) {
				library=BlibToLibraryConverter.getFile(new File(arguments.get("-l")));
			} else {
				library=null;
			}
			
			ArrayList<SearchJobData> pecanJobs=new ArrayList<SearchJobData>();
			if (diaFile.isDirectory()) {
				File[] files=diaFile.listFiles(StripeFileGenerator.getFilenameFilter());
				if (files.length==0) {
					Logger.errorLine("Your specified input (-i) directory didn't contain any .RAW files!");
					System.exit(1);
				}
				
				for (File file : files) {
					XCorDIAJobData job=new XCorDIAJobData(Optional.ofNullable(targets), Optional.ofNullable(library), file, fastaFile, factory);
					pecanJobs.add(job);
				}
			} else {
				XCorDIAJobData job=new XCorDIAJobData(Optional.ofNullable(targets), Optional.ofNullable(library), diaFile, fastaFile, factory);
				pecanJobs.add(job);
			}
			Logger.logLine("Attempting to process "+pecanJobs.size()+" searches...");
			convert(new EmptyProgressIndicator(), pecanJobs, outputFile, outputFormat, alignBetweenFiles);
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}

	public static void convertPecan(HashMap<String, String> arguments) {
		if (!arguments.containsKey("-i")||!arguments.containsKey("-f")||!arguments.containsKey("-o")) {
			Logger.errorLine("You are required to specify an input file or directory (-i), an input fasta file (-f) and an output library file (-o)");
			System.exit(1);
		}

		File diaFile=new File(arguments.get("-i"));
		File fastaFile=new File(arguments.get("-f"));
		File outputFile=new File(arguments.get("-o"));
		boolean alignBetweenFiles=ParsingUtils.getBoolean("-a", arguments, true);
		boolean writeBlib=ParsingUtils.getBoolean("-blib", arguments, false);

		final OutputFormat outputFormat = writeBlib ? OutputFormat.BLIB : OutputFormat.ELIB;

		PecanSearchParameters parameters=PecanParameterParser.parseParameters(arguments);
		PecanScoringFactory factory=new PecanOneScoringFactory(parameters, outputFile);
		Logger.logLine("SearchToLIB Pecan version "+ProgramType.getGlobalVersion().toString());

		Logger.timelessLogLine("Parameters:");
		Logger.timelessLogLine(" -i "+diaFile.getAbsolutePath());
		Logger.timelessLogLine(" -f "+fastaFile.getAbsolutePath());
		Logger.timelessLogLine(" -o "+outputFile.getAbsolutePath());
		Logger.timelessLogLine(" -a "+alignBetweenFiles);
		Logger.timelessLogLine(" -blib "+writeBlib);
		Logger.timelessLogLine(parameters.toString());

		try {
			ArrayList<FastaPeptideEntry> targets;
			if (arguments.containsKey(Pecanpie.TARGET_FASTA_TAG)) {
				targets=FastaReader.readPeptideFasta(new File(arguments.get(Pecanpie.TARGET_FASTA_TAG)), parameters);
			} else {
				targets=null;
			}
			
			ArrayList<SearchJobData> pecanJobs=new ArrayList<SearchJobData>();
			if (diaFile.isDirectory()) {
				File[] files=diaFile.listFiles(StripeFileGenerator.getFilenameFilter());
				if (files.length==0) {
					Logger.errorLine("Your specified input (-i) directory didn't contain any .RAW files!");
					System.exit(1);
				}
				
				for (File file : files) {
					PecanJobData job=new PecanJobData(Optional.ofNullable(targets), file, fastaFile, factory);
					pecanJobs.add(job);
				}
			} else {
				PecanJobData job=new PecanJobData(Optional.ofNullable(targets), diaFile, fastaFile, factory);
				pecanJobs.add(job);
			}
			Logger.logLine("Attempting to process "+pecanJobs.size()+" searches...");
			convert(new EmptyProgressIndicator(), pecanJobs, outputFile, outputFormat, alignBetweenFiles);
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}

	public static void convertEncyclopedia(HashMap<String, String> arguments) {
		if (!arguments.containsKey("-i")||!arguments.containsKey("-l")||!arguments.containsKey("-o")||!arguments.containsKey("-f")) {
			Logger.errorLine("You are required to specify an input file or directory (-i), an input library file (-l), a fasta database (-f), and an output library file (-o)");
			System.exit(1);
		}

		File diaFile=new File(arguments.get("-i"));
		File fastaFile=new File(arguments.get("-f"));
		File libraryFile=new File(arguments.get("-l"));
		File outputFile=new File(arguments.get("-o"));

		final boolean alignBetweenFiles=ParsingUtils.getBoolean("-a", arguments, true);
		final boolean writeBlib=ParsingUtils.getBoolean("-blib", arguments, false);
		final boolean alignOnly = ParsingUtils.getBoolean("-alignOnly", arguments, false);

		final SearchParameters parameters=SearchParameterParser.parseParameters(arguments);

		final OutputFormat outputFormat;

		if (alignOnly) {
			if (!alignBetweenFiles) {
				Logger.errorLine("-alignOnly requires alignment to be enabled; try running with `-a true`");
				System.exit(1);
			}

			if (writeBlib) {
				Logger.errorLine("-alignOnly requires ELIB output; try running with `-blib false`");
				System.exit(1);
			}

			if (!parameters.isQuantifySameFragmentsAcrossSamples()) {
				Logger.errorLine("-alignOnly requires -quantifyAcrossSamples true");
				System.exit(1);
			}

			outputFormat = OutputFormat.ALIB;
		} else {
			outputFormat = writeBlib ? OutputFormat.BLIB : OutputFormat.ELIB;
		}

		LibraryScoringFactory factory=new EncyclopediaOneScoringFactory(parameters);
		Logger.timelessLogLine("SearchToLIB EncyclopeDIA version "+ProgramType.getGlobalVersion().toString());

		Logger.timelessLogLine("Parameters:");
		Logger.timelessLogLine(" -i "+diaFile.getAbsolutePath());
		Logger.timelessLogLine(" -f "+fastaFile.getAbsolutePath());
		Logger.timelessLogLine(" -l "+libraryFile.getAbsolutePath());
		Logger.timelessLogLine(" -o "+outputFile.getAbsolutePath());
		Logger.timelessLogLine(" -a "+alignBetweenFiles);
		Logger.timelessLogLine(" -blib "+writeBlib);
		Logger.timelessLogLine(" -alignOnly " + alignOnly);
		Logger.timelessLogLine(parameters.toString());

		try {
			LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
			
			ArrayList<SearchJobData> pecanJobs=new ArrayList<SearchJobData>();
			if (diaFile.isDirectory()) {
				File[] files=diaFile.listFiles(StripeFileGenerator.getFilenameFilter());
				
				if (files.length==0) {
					Logger.errorLine("Your specified input (-i) directory didn't contain any .RAW files!");
					System.exit(1);
				}
				for (File file : files) {
					EncyclopediaJobData job=new EncyclopediaJobData(file, fastaFile, library, factory);
					pecanJobs.add(job);
				}
			} else {
				EncyclopediaJobData job=new EncyclopediaJobData(diaFile, fastaFile, library, factory);
				pecanJobs.add(job);
			}
			Logger.logLine("Attempting to process "+pecanJobs.size()+" searches...");

			if (arguments.containsKey("-alignmentFrom")) {
				//TODO: compute passing peptides, inferrer
				final Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides = null;
				final PeakLocationInferrerInterface inferrer = null;

				convertElibQuantOnly(new EmptyProgressIndicator(), pecanJobs, outputFile, passingPeptides, inferrer, parameters);
			} else {
				convert(new EmptyProgressIndicator(), pecanJobs, outputFile, outputFormat, alignBetweenFiles);
			}
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}

	public enum OutputFormat {
		/**
		 * Write to the ELIB format. If {@code inferrer} is present, the resulting file will be a "quantitative" ELIB,
		 * using the precomputed top-N transitions for quantification and inferred (aligned) RTs when the peptide was
		 * not detected in the initial single-file search. Additionally, quantitative matrices for peptides and proteins
		 * will be written.
		 *
		 * All passing peptides will be included.
		 *
		 * {@code globalPercolatorPeptides} should be provided if converting more than a single search, but should otherwise be empty.
		 */
		ELIB {
			@Override
			void convert(ProgressIndicator progress, List<? extends SearchJobData> jobs, File outputFile, Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides, Optional<PercolatorExecutionData> globalPercolatorFiles, Optional<PeakLocationInferrerInterface> inferrer, SearchParameters parameters) {
				convertElib(progress, jobs, outputFile, Optional.of(passingPeptides), globalPercolatorFiles, inferrer, parameters);
			}
		},

		/**
		 * Write results to the "alignment-only library" format, which records the passing peptides, RT alignment, and
		 * refined transitions for the experiment to a library file without performing any additional work. The resulting
		 * file can then be used to quantify the same targets in later separate runs of one (or more) sample(s).
		 *
		 * Note that {@code inferrer} must be present to support this export type.
		 */
		ALIB {
			@Override
			void convert(ProgressIndicator progress, List<? extends SearchJobData> jobs, File outputFile, Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides, Optional<PercolatorExecutionData> globalPercolatorFiles, Optional<PeakLocationInferrerInterface> inferrer, SearchParameters parameters) {
				if (!inferrer.isPresent()) {
					throw new IllegalArgumentException("Unable to export alignment-only library without RT alignment and transition refinement!");
				}

				convertAlib(progress, jobs, outputFile, passingPeptides, globalPercolatorFiles, inferrer.get(), parameters);
			}
		},

		/**
		 * Write to the BLIB format, suitable for use with Skyline. Only quantifiable peptides will be written from each
		 * search. Additionally, a TSV "integration" file will be written with details of the peptides included in the library.
		 */
		BLIB {
			@Override
			void convert(ProgressIndicator progress, List<? extends SearchJobData> jobs, File outputFile, Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides, Optional<PercolatorExecutionData> globalPercolatorFiles, Optional<PeakLocationInferrerInterface> inferrer, SearchParameters parameters) {
				convertBlib(progress, jobs, outputFile, Optional.of(passingPeptides.x), inferrer);
			}
		};

		/**
		 * Write results to the given location in this format. Typically this method should only be called from
		 * {@link SearchToBLIB#convert(ProgressIndicator, List, File, OutputFormat, boolean)} which will handle either
		 * reading or computing the necessary information for a group of samples.
		 *
		 * Will also compute and output related information in some cases, depending on the format.
		 *
		 * @param progress A progress indicator that will be used during the conversion process
		 * @param jobs The jobs whose results should be included in the output file
		 * @param outputFile The location where the new library will be created (will be overwritten if it exists)
		 * @param passingPeptides The results of running Percolator to determine the list of peptides that will be
		 *                        included, as returned by {@link PercolatorReader#getPassingPeptidesFromTSV}
		 * @param globalPercolatorFiles Used by some formats to get additional information when Percolator has been run on
		 *                              results from multiple input files
		 * @param inferrer If aligning between files, the inferrer which provides RT alignment and consistent, refined transitions
		 * @param parameters The parameters that should be used during conversion and (in some cases) written to the output file
		 */
		abstract void convert(ProgressIndicator progress, List<? extends SearchJobData> jobs, File outputFile, Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides, Optional<PercolatorExecutionData> globalPercolatorFiles, Optional<PeakLocationInferrerInterface> inferrer, SearchParameters parameters);
	}

	/**
	 * Legacy form of {@link #convert(ProgressIndicator, List, File, OutputFormat, boolean)} which supports only
	 * ELIB and BLIB formats.
	 *
	 * @see #convert(ProgressIndicator, List, File, OutputFormat, boolean)
	 *
	 * @deprecated it's better to directly specify the desired output format with an enum constant
	 */
	@Deprecated
	public static void convert(ProgressIndicator progress, List<? extends SearchJobData> pecanJobs, File libFile, boolean writeBlib, boolean alignBetweenFiles) {
		convert(
				progress,
				pecanJobs,
				libFile,
				writeBlib ? OutputFormat.BLIB : OutputFormat.ELIB,
				alignBetweenFiles
		);
	}

	/**
	 * For the given previously-run single-file searches (jobs), gather or compute the necessary information to create
	 * a combined output in the given format. This handles the core jobs of (if necessary) running Percolator, reading
	 * Percolator results to determine the set of global passing peptides, performing retention time alignment and
	 * transition refinement (if {@code alignBetweenFiles} is true), and writing results to the given output file, which
	 * may involve additional work like quantifying peptides in each single file.
	 *
	 * @param progress A progress indicator that will be used during the conversion process
	 * @param pecanJobs The jobs whose results should be included in the output file
	 * @param libFile The location where the new library will be created (will be overwritten if it exists)
	 * @param outputFormat The format which should be written
	 * @param alignBetweenFiles If RT alignment
	 */
	public static void convert(ProgressIndicator progress, List<? extends SearchJobData> pecanJobs, File libFile, OutputFormat outputFormat, boolean alignBetweenFiles) {
		ArrayList<SearchJobData> processedJobs=new ArrayList<SearchJobData>();
		ArrayList<File> featureFiles=new ArrayList<File>();
		SearchJobData representativeJob=null;

		pecanJobs = Lists.newArrayList(pecanJobs); // mutable copy

		// Sort files in alphabetical order for deterministic Percolator sampling
		Collections.sort(pecanJobs, (a, b) -> a.getDiaFileReader().getOriginalFileName().compareTo(b.getDiaFileReader().getOriginalFileName()));
		
		for (int i=0; i<pecanJobs.size(); i++) {
			SearchJobData job=pecanJobs.get(i);
			if (!job.hasBeenRun()) {
				Logger.logLine("Can't find a "+job.getSearchType()+" analysis of "+job.getDiaFileReader().getOriginalFileName()+", skipping extraction on that file.");
				continue;
			} else {
				processedJobs.add(job);
			}
			if (representativeJob==null) {
				representativeJob=job;
			}
			featureFiles.add(job.getPercolatorFiles().getInputTSV());
		}
		pecanJobs=processedJobs;

		if (representativeJob==null) {
			Logger.errorLine("Can't find any representative jobs! Failing...");

			for (int i=0; i<pecanJobs.size(); i++) {
				SearchJobData job=pecanJobs.get(i);
				Logger.errorLine(" Checking raw file "+(i+1)+": "+job.getDiaFileReader().getFile().exists());
				Logger.errorLine(" Checking feature file "+(i+1)+": "+job.getPercolatorFiles().getInputTSV().exists());
				Logger.errorLine(" Checking result file "+(i+1)+": "+job.getPercolatorFiles().getPeptideOutputFile().exists());
			}
			return;
		}
		Logger.logLine("Using "+representativeJob.getDiaFileReader().getOriginalFileName()+" to extract representative search parameters");
		SearchParameters parameters=representativeJob.getParameters();

		String filename=libFile.getName();
		if (filename.lastIndexOf('.')>0) {
			filename=filename.substring(0, filename.lastIndexOf('.'));
		}
		File bigFeatureFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_features.txt");
		File bigPercolatorFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_results.txt");
		File bigPercolatorDecoyFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_decoy.txt");
		File bigPercolatorProteinFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_protein_results.txt");
		File bigPercolatorProteinDecoyFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_protein_decoy.txt");
		PercolatorExecutionData bigPercolatorFiles=new PercolatorExecutionData(bigFeatureFile, representativeJob.getPercolatorFiles().getFastaFile(), bigPercolatorFile, bigPercolatorDecoyFile, bigPercolatorProteinFile, bigPercolatorProteinDecoyFile, parameters);
		
		final float threshold=parameters.getEffectivePercolatorThreshold();
		try {
			Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides;
			boolean runningPercolator=true;
			if (featureFiles.size()==1) {
				Logger.logLine("Only one file, so no need to re-run Percolator.");
				// if there's only one file then don't need to re-run percolator
				passingPeptides=PercolatorReader.getPassingPeptidesFromTSV(representativeJob.getPercolatorFiles().getPeptideOutputFile(), parameters, false);
				runningPercolator=false;
				
			} else if (parameters.isDoNotUseGlobalFDR()) {
				Logger.logLine("Warning, user asked to not use global FDR!");
				passingPeptides=getPeptidesWithoutGlobalFDR(pecanJobs, parameters).x;
				runningPercolator=false;
				
			} else if (bigPercolatorFile.exists()&&bigPercolatorFile.canRead()&&bigPercolatorDecoyFile.exists()&&bigPercolatorDecoyFile.canRead()) {
				Logger.logLine("Found previously run global Percolator.");
				// if we've already run percolator then don't need to re-run percolator
				passingPeptides=PercolatorReader.getPassingPeptidesFromTSV(bigPercolatorFile, parameters, false);
				runningPercolator=false;
			} else {
				Logger.logLine("Running global Percolator analysis.");
				TableConcatenator.concatenateTables(featureFiles, bigFeatureFile);
				
				// delete if exists
				if (bigPercolatorFiles.getModelFile().exists()) {
					bigPercolatorFiles.getModelFile().delete();
				}
				int modelNumber = Integer.MAX_VALUE; // always use the last model (if reusing a model)
				passingPeptides=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorVersionNumber(), bigPercolatorFiles, threshold, parameters.getAAConstants(), modelNumber);
			}

			Logger.logLine("Identified "+passingPeptides.x.size()+" peptides across all files at a "+(threshold*100.0f)+"% FDR threshold.");

			boolean foundLibrary=false;
			if ((!runningPercolator)&&libFile.exists()&&libFile.canRead()) {
				// didn't have to run percolator, so check if we can read the lib file
				final LibraryFile lib=new LibraryFile();
				try {
					lib.openFile(libFile);
					Logger.logLine("Found library file and tested for reading. It seems ok so proceeding with that file!");
					foundLibrary=true;
				} catch (Exception e) {
					Logger.logLine("Found library file and tested for reading. Reading failed, so overwriting!");
				} finally {
					lib.close();
				}
			}
			
			if (!foundLibrary) {
				Optional<PeakLocationInferrerInterface> inferrer;
				if (alignBetweenFiles) {
					Logger.logLine("Inferring peak boundaries across files...");
					try {
						inferrer=Optional.of(AlternatePeakLocationInferrer.getAlignmentData(new EmptyProgressIndicator(), pecanJobs, passingPeptides.x, parameters));
						Logger.logLine("...Finished peak inference.");
					} catch (Exception e) {
						Logger.errorLine("RT alignment between files failed! Perhaps this is to build a chromatogram library and not a quantitative experiment? Attempting to recover without alignment.");
						inferrer=Optional.empty();
					}
				} else {
					Logger.logLine("No RT alignment between files necessary.");
					inferrer=Optional.empty();
				}

				outputFormat.convert(progress, pecanJobs, libFile, passingPeptides, Optional.ofNullable(featureFiles.size() == 1 ? null : bigPercolatorFiles), inferrer, parameters);
			}
			progress.update(passingPeptides.x.size()+" peptides identified at "+(threshold*100.0f)+"% FDR", 1.0f);
		} catch (IOException ioe) {
			Logger.errorLine("Error creating concatenated feature file");
			Logger.errorException(ioe);
		} catch (InterruptedException ie) {
			Logger.errorLine("Error creating concatenated feature file");
			Logger.errorException(ie);
		}
	}

	private static Pair<Pair<ArrayList<PercolatorPeptide>, Float>, Pair<ArrayList<PercolatorPeptide>, Float>> getPeptidesWithoutGlobalFDR(List<? extends SearchJobData> pecanJobs, SearchParameters parameters) {
		Pair<ArrayList<PercolatorPeptide>, Float> resultDecoyPeptides=new Pair<ArrayList<PercolatorPeptide>, Float>(new ArrayList<>(), -1f);
		HashMap<String, ScoredObject<PeptidePrecursor>> decoyMap=new HashMap<>();
		for (SearchJobData job : pecanJobs) {
			ArrayList<PercolatorPeptide> individualSamplePeptides=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideDecoyFile(), parameters, false).x;
			for (PercolatorPeptide peptide : individualSamplePeptides) {
				ScoredObject<PeptidePrecursor> obj=decoyMap.get(peptide.getPeptideModSeq());
				if (obj==null||obj.x>peptide.getPosteriorErrorProb()) {
					decoyMap.put(peptide.getPeptideModSeq(), new ScoredObject<PeptidePrecursor>(peptide.getPosteriorErrorProb(), peptide));
				}
			}
		}
		for (ScoredObject<PeptidePrecursor> precursor : decoyMap.values()) {
			resultDecoyPeptides.x.add((PercolatorPeptide)precursor.y);
		}
		
		Pair<ArrayList<PercolatorPeptide>, Float> resultTargetPeptides=new Pair<ArrayList<PercolatorPeptide>, Float>(new ArrayList<>(), -1f);
		HashMap<String, ScoredObject<PeptidePrecursor>> targetMap=new HashMap<>();
		for (SearchJobData job : pecanJobs) {
			ArrayList<PercolatorPeptide> individualSamplePeptides=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideOutputFile(), parameters, false).x;
			for (PercolatorPeptide peptide : individualSamplePeptides) {
				ScoredObject<PeptidePrecursor> obj=targetMap.get(peptide.getPeptideModSeq());
				if (obj==null||obj.x>peptide.getPosteriorErrorProb()) {
					targetMap.put(peptide.getPeptideModSeq(), new ScoredObject<PeptidePrecursor>(peptide.getPosteriorErrorProb(), peptide));
				}
			}
		}
		for (ScoredObject<PeptidePrecursor> precursor : targetMap.values()) {
			resultTargetPeptides.x.add((PercolatorPeptide)precursor.y);
		}
		
		return new Pair<Pair<ArrayList<PercolatorPeptide>,Float>, Pair<ArrayList<PercolatorPeptide>,Float>>(resultTargetPeptides, resultDecoyPeptides);
	}
	
	static void convertBlib(ProgressIndicator progress, List<? extends SearchJobData> pecanJobs, File blibFile, Optional<ArrayList<PercolatorPeptide>> passingPeptides, Optional<PeakLocationInferrerInterface> inferrer) {
		try {
			BlibFile blib=new BlibFile();
			blib.openFile();
			blib.setUserFile(blibFile);
			blib.dropIndices();
			int[] counterTotals=new int[] {0,0,0};
			
			File integrationFile=new File(blibFile.getAbsolutePath()+".integration.txt");
			PrintWriter integrationFileWriter=new PrintWriter(integrationFile, "UTF-8");
			integrationFileWriter.println("File Name\tPeptide Modified Sequence\tMin Start Time\tMax End Time\tPrecursor Charge\tPrecursorIsDecoy\tIon Count\tRetention Time Center\tTIC");

			float increment=1.0f/pecanJobs.size();
			for (int i=0; i<pecanJobs.size(); i++) {
				SearchJobData job=pecanJobs.get(i);
				if (!job.hasBeenRun()) {
					continue;
				}
				ProgressIndicator subProgress=new SubProgressIndicator(progress, increment);
				
				ArrayList<PercolatorPeptide> globalPassingPeptides;
				ArrayList<PercolatorPeptide> localPassingPeptides=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideOutputFile(), pecanJobs.get(i).getParameters(), false).x;
				if (passingPeptides.isPresent()) {
					globalPassingPeptides=passingPeptides.get();
				} else {
					globalPassingPeptides=localPassingPeptides;
				}
				
				counterTotals=convertFileBlib(subProgress, job, globalPassingPeptides, localPassingPeptides, counterTotals, inferrer, integrationFileWriter, blib);
			}
			integrationFileWriter.flush();
			integrationFileWriter.close();

			blib.createIndices();
			blib.saveFile();
			blib.close();
		} catch (IOException ioe) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(ioe);
		} catch (SQLException sqle) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(sqle);
		}
	}

	/**
	 * trims to quantifiable peptides! for loading into skyline!
	 */
	static int[] convertFileBlib(ProgressIndicator subProgress, SearchJobData job, ArrayList<PercolatorPeptide> globalPassingPeptides, ArrayList<PercolatorPeptide> localPassingPeptides, int[] counterTotals, Optional<PeakLocationInferrerInterface> inferrer, PrintWriter integrationFileWriter, BlibFile blib) throws IOException, SQLException {
		final String diaFileName = job.getDiaFileReader().getOriginalFileName();

		Logger.logLine("Reading Percolator Results from "+ diaFileName +"...");
		subProgress.update(diaFileName +": Reading Percolator Results", 0.0f);

		final StripeFileInterface stripeFile = job.getDiaFileReader();

		Logger.logLine("Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides from "+ diaFileName +"...");
		subProgress.update(diaFileName +": Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides", 0.00001f);

		LibraryInterface library=null;
		if (job instanceof EncyclopediaJobData) {
			library=((EncyclopediaJobData)job).getLibrary();
		}
		//ArrayList<IntegratedLibraryEntry> libraryEntries=SearchFeatureReader.parseSearchFeatures(featureFile, globalPassingPeptides, localPassingPeptides, stripeFile, Optional.ofNullable((LibraryFile)null), job.getParameters());
		ArrayList<IntegratedLibraryEntry> libraryEntries=PeptideQuantExtractor.parseSearchFeatures(subProgress, job, true, globalPassingPeptides, localPassingPeptides, inferrer, stripeFile, library, job.getParameters());
		stripeFile.close();
		
		for (IntegratedLibraryEntry entry : libraryEntries) {
			String peptideModSeq=PeptideUtils.formatForSkylinePeakBoundaries(entry.getPeptideModSeq());
			integrationFileWriter.println(diaFileName +"\t"+peptideModSeq+"\t"+entry.getRtRange().getStart()/60f+"\t"+entry.getRtRange().getStop()/60f+"\t"+entry.getPrecursorCharge()+"\tFALSE\t"+entry.getIonCount()+"\t"+entry.getRetentionTime()/60f+"\t"+entry.getTIC());
		}
		integrationFileWriter.flush();
		
		ArrayList<LibraryEntry> recasted=new ArrayList<LibraryEntry>();
		for (IntegratedLibraryEntry entry : libraryEntries) {
			recasted.add(entry);
		}

		Logger.logLine("Writing Skyline BLIB from "+ diaFileName +"...");
		subProgress.update(diaFileName +": Writing Skyline BLIB", 0.99999f);

		counterTotals=blib.addLibrary(job, recasted, counterTotals[0], counterTotals[1], counterTotals[2]);
		subProgress.update(diaFileName +": Finished writing to Skyline BLIB at"+new Date().toString(), 1.0f);
		return counterTotals;
	}

	static void convertElib(ProgressIndicator progress, SearchJobData pecanJob, File elibFile, SearchParameters parameters) {
		ArrayList<SearchJobData> jobs=new ArrayList<>();
		jobs.add(pecanJob);

		convertElib(progress, jobs, elibFile, Optional.empty(), Optional.empty(), Optional.empty(), parameters);
	}

	static void convertElib(ProgressIndicator progress, List<? extends SearchJobData> pecanJobs, File elibFile, Optional<Pair<ArrayList<PercolatorPeptide>, Float>> passingPeptides, Optional<PercolatorExecutionData> globalPercolatorFiles, Optional<PeakLocationInferrerInterface> inferrer, SearchParameters parameters) {
		try {
			LibraryFile elib=new LibraryFile();
			elib.openFile();
			elib.dropIndices();

			float increment=1.0f/pecanJobs.size();
			for (int i=0; i<pecanJobs.size(); i++) {
				SearchJobData job=pecanJobs.get(i);
				if (!job.hasBeenRun()) {
					continue;
				}
				ProgressIndicator subProgress=new SubProgressIndicator(progress, increment);
				
				ArrayList<PercolatorPeptide> globalPassingPeptides;
				Pair<ArrayList<PercolatorPeptide>, Float> localPassingPeptides=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideOutputFile(), pecanJobs.get(i).getParameters(), false);
				if (passingPeptides.isPresent()) {
					globalPassingPeptides=passingPeptides.get().x;
				} else {
					globalPassingPeptides=localPassingPeptides.x;
				}

				Logger.logLine(job.getDiaFileReader().getOriginalFileName()+": Number of global peptides: "+globalPassingPeptides.size()+" vs local peptides: "+localPassingPeptides.x.size());
				
				convertFileElib(subProgress, job, globalPassingPeptides, localPassingPeptides.x, inferrer, elib, pecanJobs.size()>1);

				if ((!globalPercolatorFiles.isPresent())) {
					if (job.hasBeenRun()) {
						Pair<ArrayList<PercolatorPeptide>, Float> targets=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideOutputFile(), parameters, true);
						Pair<ArrayList<PercolatorPeptide>, Float> decoys=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideDecoyFile(), parameters, true);
						Logger.logLine("Writing local target/decoy peptides: "+targets.x.size()+"/"+decoys.x.size()+", pi0: "+targets.y);
						elib.addTargetDecoyPeptides(targets.x, decoys.x);
						elib.addMetadata("pi0", Float.toString(targets.y));
						elib.addProteinsFromPercolator(targets.x);
						elib.addProteinsFromPercolator(decoys.x);
						
						Pair<ArrayList<PercolatorProteinGroup>, ArrayList<PercolatorProteinGroup>> targetDecoyProteins=ParsimonyProteinGrouper.groupProteins(targets.x, decoys.x, parameters.getPercolatorProteinThreshold(), parameters.getAAConstants());
						Logger.logLine("Writing local target/decoy proteins: "+targetDecoyProteins.x.size()+"/"+targetDecoyProteins.y.size());
						elib.addTargetDecoyProteins(job.getDiaFileReader().getOriginalFileName(), targetDecoyProteins.x, targetDecoyProteins.y);

						job.getPercolatorFiles()
								.getPercolatorExecutableVersion()
								.ifPresent((ThrowingConsumer<String>) version -> {
									elib.addMetadata(LibraryFile.PERCOLATOR_VERSION, version);
								});
					}
				}
				
				subProgress.update("Wrote "+globalPassingPeptides.size()+" peptides identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
			}

			ArrayList<PercolatorProteinGroup> proteins=null;
			if (globalPercolatorFiles.isPresent()) {
				proteins = writePercolatorToElib(elib, globalPercolatorFiles.get(), pecanJobs, parameters);
			}

			writeElibMetadata(elib, pecanJobs, parameters, inferrer.isPresent());

			elib.createIndices();
			elib.saveAsFile(elibFile);

			if (proteins!=null) {
				if (inferrer.isPresent()) {
					try {
						ArrayList<ProteinGroupInterface> proteinGroups=new ArrayList<>();
						for (ProteinGroupInterface pg : proteins) {
							proteinGroups.add(pg);
						}
						LibraryReportExtractor.extractMatrix(elib, proteinGroups, true);
					} catch (DataFormatException e) {
						Logger.errorException(e);
					}
				}
			}

			elib.close();
		} catch (IOException ioe) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(ioe);
		} catch (SQLException sqle) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(sqle);
		}
	}

	/**
	 * Does not limit to quantifiable! Reports all potential peaks!
	 */
	static void convertFileElib(ProgressIndicator subProgress, SearchJobData job, ArrayList<PercolatorPeptide> globalPassingPeptides, ArrayList<PercolatorPeptide> localPassingPeptides, Optional<PeakLocationInferrerInterface> inferrer, LibraryFile elib, boolean combineJobs) throws IOException, SQLException {
		String diaFileName=job.getDiaFileReader().getOriginalFileName();
		Logger.logLine("Reading Percolator Results from "+diaFileName+"...");
		subProgress.update(diaFileName+": Reading Percolator Results", 0.0f);

		final StripeFileInterface stripeFile = job.getDiaFileReader();

		Logger.logLine("Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides from "+diaFileName+"...");
		subProgress.update(diaFileName+": Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides", 0.00001f);

		elib.addTIC(stripeFile);

		inferrer.ifPresent(inf -> elib.addRtAlignment(job, inf));

		ArrayList<IntegratedLibraryEntry> libraryEntries;
		if (job instanceof QuantitativeSearchJobData) {
			LibraryInterface library=null;
			if (job instanceof EncyclopediaJobData) {
				library=((EncyclopediaJobData)job).getLibrary();
			}
			libraryEntries=PeptideQuantExtractor.parseSearchFeatures(subProgress, job, false, globalPassingPeptides, localPassingPeptides, inferrer, stripeFile, library, job.getParameters());
		} else {
			HashMap<String, PSMData> targetPSMs=PeptideQuantExtractor.findTargetPSMData(job, globalPassingPeptides, localPassingPeptides, inferrer, job.getParameters());
			libraryEntries=DDAPrecursorIntegrator.integrateSearch(subProgress, targetPSMs, stripeFile, job.getParameters());
		}
		stripeFile.close();
		
		Logger.logLine("Writing Encyclopedia ELIB from "+diaFileName+" ("+libraryEntries.size()+" entries)...");
		subProgress.update(diaFileName+": Writing Encyclopedia ELIB", 0.99999f);
		
		Optional<HashMap<String, ModificationLocalizationData>> localizationData;
		if (!combineJobs&&job instanceof ThesaurusJobData) {
			Logger.logLine("Reading localization data from disk...");
			localizationData=Optional.of(LocalizationDataToTSVConsumer.readLocalizationFile(((ThesaurusJobData)job).getLocalizationFile(), globalPassingPeptides, job.getParameters()));
		} else if (!combineJobs&&job instanceof VariantXCorDIAJobData) {
			Logger.logLine("Reading localization data from disk...");
			localizationData=Optional.of(LocalizationDataToTSVConsumer.readLocalizationFile(((VariantXCorDIAJobData)job).getLocalizationFile(), globalPassingPeptides, job.getParameters()));
		} else {
			localizationData=Optional.empty();
		}

		elib.addIntegratedEntries(libraryEntries, inferrer, localizationData, job.getParameters().getAAConstants(), job.getParameters().getPercolatorThreshold());
		

		Logger.logLine("Finished writing to Encyclopedia ELIB at "+new Date().toString());
		subProgress.update(diaFileName+": Finished writing to Encyclopedia ELIB at "+new Date().toString(), 1.0f);
	}

	/**
	 * Read the set of passing peptides from Percolator and write them to the ELIB with associated metadata.
	 * Perform protein inference and write protein scores/q-values/PEPs to the ELIB.
	 *
	 * @param elib The (open) ELIB where results will be written.
	 * @param percolatorExecutionData Used to read the list of passing peptides and associated scores/metadata. If the
	 *                                global results file doesn't exist this will be ignored and the calculation will
	 *                                fall back to use {@code jobs}.
	 * @param jobs Ignored, unless global Percolator results don't exist, in which case the passing peptides are read
	 *             directly from these jobs, without global FDR control.
	 * @return The inferred set of protein groups.
	 */
	private static ArrayList<PercolatorProteinGroup> writePercolatorToElib(LibraryFile elib, PercolatorExecutionData percolatorExecutionData, List<? extends SearchJobData> jobs, SearchParameters parameters) throws IOException, SQLException {
		return writePercolatorToElib(elib, percolatorExecutionData, Optional.of(jobs), parameters);
	}

	/**
	 * Read the set of passing peptides from Percolator and write them to the ELIB with associated metadata.
	 * Perform protein inference and write protein scores/q-values/PEPs to the ELIB.
	 *
	 * @param elib The (open) ELIB where results will be written.
	 * @param percolatorExecutionData Used to read the list of passing peptides and associated scores/metadata.
	 * @return The inferred set of protein groups.
	 */
	private static ArrayList<PercolatorProteinGroup> writePercolatorToElib(LibraryFile elib, PercolatorExecutionData percolatorExecutionData, SearchParameters parameters) throws IOException, SQLException {
		return writePercolatorToElib(elib, percolatorExecutionData, Optional.empty(), parameters);
	}

	/**
	 * Read the set of passing peptides from Percolator and write them to the ELIB with associated metadata.
	 * Perform protein inference and write protein scores/q-values/PEPs to the ELIB.
	 *
	 * @param elib The (open) ELIB where results will be written.
	 * @param percolatorExecutionData Used to read the list of passing peptides and associated scores/metadata. If the
	 *                                global results file doesn't exist this will be ignored and the calculation will
	 *                                fall back to use {@code jobs}.
	 * @param jobs Ignored, unless global Percolator results don't exist, in which case the passing peptides are read
	 *             directly from these jobs, without global FDR control. An exception will be raised if this fallback
	 *             is necessary but {@code jobs} is not present.
	 * @return The inferred set of protein groups.
	 */
	private static ArrayList<PercolatorProteinGroup> writePercolatorToElib(LibraryFile elib, PercolatorExecutionData percolatorExecutionData, Optional<List<? extends SearchJobData>> jobs, SearchParameters parameters) throws IOException, SQLException {
		Pair<ArrayList<PercolatorPeptide>, Float> targets=null;
		Pair<ArrayList<PercolatorPeptide>, Float> decoys=null;
		if (percolatorExecutionData.getPeptideOutputFile().exists()) {
			targets=PercolatorReader.getPassingPeptidesFromTSV(percolatorExecutionData.getPeptideOutputFile(), parameters, true);
			decoys=PercolatorReader.getPassingPeptidesFromTSV(percolatorExecutionData.getPeptideDecoyFile(), parameters, true);
		} else if (jobs.isPresent()) {
			Pair<Pair<ArrayList<PercolatorPeptide>, Float>, Pair<ArrayList<PercolatorPeptide>, Float>> withoutFDR=getPeptidesWithoutGlobalFDR(jobs.get(), parameters);
			targets=withoutFDR.x;
			decoys=withoutFDR.y;
		} else {
			throw new IllegalStateException("Unable to get passing peptides: no global Percolator results file or individual jobs!");
		}

		Logger.logLine("Writing global target/decoy peptides: "+targets.x.size()+"/"+decoys.x.size()+", pi0: "+targets.y);
		elib.addTargetDecoyPeptides(targets.x, decoys.x);
		elib.addMetadata("pi0", Float.toString(targets.y));
		elib.addProteinsFromPercolator(targets.x);
		elib.addProteinsFromPercolator(decoys.x);

		Pair<ArrayList<PercolatorProteinGroup>, ArrayList<PercolatorProteinGroup>> targetDecoyProteins=ParsimonyProteinGrouper.groupProteins(targets.x, decoys.x, parameters.getPercolatorProteinThreshold(), parameters.getAAConstants());

		Logger.logLine("Writing global target/decoy proteins: "+targetDecoyProteins.x.size()+"/"+targetDecoyProteins.y.size());
		elib.addTargetDecoyProteins("global", targetDecoyProteins.x, targetDecoyProteins.y);

		percolatorExecutionData
				.getPercolatorExecutableVersion()
				.ifPresent((ThrowingConsumer<String>) version -> {
					elib.addMetadata(LibraryFile.PERCOLATOR_VERSION, version);
				});

		return targetDecoyProteins.x;
	}

	private static void writeElibMetadata(LibraryFile elib, List<? extends SearchJobData> jobs, SearchParameters parameters, boolean align) throws IOException, SQLException {
		final HashMap<String, String> parameterMap = parameters.toParameterMap();
		parameterMap.put("RT align between samples", Boolean.toString(align));
		for (int i = 0; i < jobs.size(); i++) {
			final SearchJobData job = jobs.get(i);
			parameterMap.put(job.getDiaFileReader().getOriginalFileName() + " search type", job.getSearchType());
			if (job instanceof EncyclopediaJobData) {
				parameterMap.put(job.getDiaFileReader().getOriginalFileName() + " library", ((EncyclopediaJobData) job).getLibrary().getName());
			} else if (job instanceof PecanJobData) {
				parameterMap.put(job.getDiaFileReader().getOriginalFileName() + " fasta", ((PecanJobData) job).getFastaFile().getName());
				parameterMap.put(job.getDiaFileReader().getOriginalFileName() + " used narrow target list", Boolean.toString(((PecanJobData) job).getTargetList().isPresent()));
			} else if (job instanceof XCorDIAJobData) {
				Optional<LibraryInterface> maybeLibrary = ((XCorDIAJobData) job).getLibrary();
				if (maybeLibrary.isPresent()) {
					parameterMap.put(job.getDiaFileReader().getOriginalFileName() + " library", maybeLibrary.get().getName());
				}
				parameterMap.put(job.getDiaFileReader().getOriginalFileName() + " fasta", ((XCorDIAJobData) job).getFastaFile().getName());
				parameterMap.put(job.getDiaFileReader().getOriginalFileName() + " used narrow target list", Boolean.toString(((XCorDIAJobData) job).getTargetList().isPresent()));
			}
		}
		elib.addMetadata(parameterMap);

		elib.setSources(jobs);
	}

	private static void convertAlib(ProgressIndicator progress, List<? extends SearchJobData> jobs, File outputFile, Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides, Optional<PercolatorExecutionData> globalPercolatorFiles, PeakLocationInferrerInterface inferrer, SearchParameters parameters) {
		if (Objects.requireNonNull(jobs, "No jobs provided").isEmpty()) {
			throw new IllegalArgumentException("No jobs provided");
		}

		if (!parameters.isQuantifySameFragmentsAcrossSamples()) {
			throw new IllegalArgumentException("Unable to export alignment-only library without -quantifyAcrossSamples!");
		}

		if (!globalPercolatorFiles.isPresent()) {
			if (jobs.size() == 1) {
				globalPercolatorFiles = Optional.of(jobs.iterator().next().getPercolatorFiles());
			} else {
				throw new IllegalArgumentException("Global percolator files must be provided for more than one job!");
			}
		}

		try {
			final LibraryFile elib = new LibraryFile();
			try {
				elib.openFile();
				elib.dropIndices();

				float increment = 1.0f / jobs.size();
				for (int i = 0; i < jobs.size(); i++) {
					final SearchJobData job = jobs.get(i);

					final ProgressIndicator subProgress = new SubProgressIndicator(progress, increment);

					if (!job.hasBeenRun()) {
						subProgress.update("Skipping incomplete job: " + job.getDiaFileReader().getOriginalFileName(), 1f);
						continue;
					}

					elib.addRtAlignment(job, inferrer);

//					elib.addEntries(job.getR); //TODO: write entries for passing peptides from this job

//					subProgress.update("Wrote "+passingPeptides.x.size()+" peptides identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
				}

				final PercolatorExecutionData percolatorExecutionData = globalPercolatorFiles.get();
				if (!percolatorExecutionData.getPeptideOutputFile().exists()) {
					throw new IllegalArgumentException("Could not read Percolator results!", new FileNotFoundException(percolatorExecutionData.getPeptideOutputFile().getAbsolutePath()));
				}

				writePercolatorToElib(elib, percolatorExecutionData, parameters);

				writeElibMetadata(elib, jobs, parameters, true); // align is required for ALIB

				elib.createIndices();
				elib.saveAsFile(outputFile);
			} finally {
				elib.close();
			}
		} catch (IOException | SQLException ioe) {
			Logger.errorLine("Error creating ELIB file");
			Logger.errorException(ioe);
			throw new EncyclopediaException("Error creating ELIB file", ioe);
		}
	}

	static void convertElibQuantOnly(ProgressIndicator progress, List<? extends SearchJobData> pecanJobs, File elibFile, Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides, PeakLocationInferrerInterface inferrer, SearchParameters parameters) {
		try {
			LibraryFile elib=new LibraryFile();
			elib.openFile();
			elib.dropIndices();

			float increment=1.0f/pecanJobs.size();
			for (int i=0; i<pecanJobs.size(); i++) {
				SearchJobData job = pecanJobs.get(i);
				if (!job.hasBeenRun()) {
					Logger.errorLine("Unable to process " + job.getDiaFileReader().getOriginalFileName() + " because its results are missing. Continuing.");
					continue;
				}
				ProgressIndicator subProgress = new SubProgressIndicator(progress, increment);

				ArrayList<PercolatorPeptide> globalPassingPeptides = passingPeptides.x;
				Pair<ArrayList<PercolatorPeptide>, Float> localPassingPeptides = PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideOutputFile(), pecanJobs.get(i).getParameters(), false);

				Logger.logLine(job.getDiaFileReader().getOriginalFileName() + ": Number of global peptides: " + globalPassingPeptides.size() + " vs local peptides: " + localPassingPeptides.x.size());

				convertFileElib(subProgress, job, globalPassingPeptides, localPassingPeptides.x, Optional.of(inferrer), elib, pecanJobs.size() > 1);
			}

			//TODO: get proteins as argument
			ArrayList<PercolatorProteinGroup> proteins=null;

			writeElibMetadata(elib, pecanJobs, parameters, true);

			elib.createIndices();
			elib.saveAsFile(elibFile);

			Objects.requireNonNull(proteins, "Unable to proceed without previously-computed protein groups!");

			try {
				ArrayList<ProteinGroupInterface> proteinGroups=new ArrayList<>();
				for (ProteinGroupInterface pg : proteins) {
					proteinGroups.add(pg);
				}
				LibraryReportExtractor.extractMatrix(elib, proteinGroups, true);
			} catch (DataFormatException e) {
				Logger.errorException(e);
			}

			elib.close();
		} catch (IOException ioe) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(ioe);
		} catch (SQLException sqle) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(sqle);
		}
	}
}
