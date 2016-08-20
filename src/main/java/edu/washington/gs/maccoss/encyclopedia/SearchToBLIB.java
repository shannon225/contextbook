package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ParsimonyProteinGrouper;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideQuantExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.LibraryReportExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableConcatenator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;

public class SearchToBLIB {
	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		if (arguments.size()==0) {
			SearchGUIMain.runGUI(false);
			
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("SearchToLIB Help");
			Logger.logLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.logLine("Required Parameters: ");
			Logger.logLine("\t-i\tinput .DIA or .MZML file or directory");
			Logger.logLine("\t-l\toriginal library .ELIB file");
			Logger.logLine("\t-o\toutput library .ELIB file");
			Logger.logLine("\t-a\talign between files (default=true)");
			
			System.exit(1);
			
		} else if (arguments.containsKey("-v")||arguments.containsKey("-version")||arguments.containsKey("--version")) {
			Logger.logLine("Encyclopedia SearchToLIB version "+EncyclopediaOneScoringFactory.version);
			System.exit(1);
			
		} else {
			if (!arguments.containsKey("-i")||!arguments.containsKey("-l")||!arguments.containsKey("-o")) {
				Logger.errorLine("You are required to specify an input file or directory (-i), an input library file (-l) and an output library file (-o)");
				System.exit(1);
			}

			File diaFile=new File(arguments.get("-i"));
			File libraryFile=new File(arguments.get("-l"));
			File outputFile=new File(arguments.get("-o"));
			boolean alignBetweenFiles=SearchParameterParser.getBoolean("-a", arguments, true);
			
			SearchParameters parameters=SearchParameterParser.parseParameters(arguments);
			LibraryScoringFactory factory=new EncyclopediaOneScoringFactory(parameters);
			Logger.logLine("Encyclopedia version "+factory.getVersion());

			Logger.logLine("Parameters:");
			Logger.logLine(" -i "+diaFile.getAbsolutePath());
			Logger.logLine(" -l "+libraryFile.getAbsolutePath());
			Logger.logLine(" -o "+outputFile.getAbsolutePath());
			Logger.logLine(" -a"+alignBetweenFiles);
			Logger.logLine(parameters.toString());

			try {
				LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
				
				ArrayList<SearchJobData> pecanJobs=new ArrayList<SearchJobData>();
				if (diaFile.isDirectory()) {
					File[] files=diaFile.listFiles(new SimpleFilenameFilter(MzmlToDIAConverter.MZML_EXTENSION));
					if (files.length==0) {
						files=diaFile.listFiles(new SimpleFilenameFilter(StripeFile.DIA_EXTENSION));
					}
					if (files.length==0) {
						Logger.errorLine("Your specified input (-i) directory didn't contain any .mzML or .DIA files!");
						System.exit(1);
					}
					for (File file : files) {
						EncyclopediaJobData job=new EncyclopediaJobData(file, library, factory);
						pecanJobs.add(job);
					}
				} else {
					EncyclopediaJobData job=new EncyclopediaJobData(diaFile, library, factory);
					pecanJobs.add(job);
				}
				convert(new EmptyProgressIndicator(), pecanJobs, outputFile, false, alignBetweenFiles);
			} catch (Exception e) {
				System.err.println("Encountered Fatal Error!");
				e.printStackTrace();
			}
		}
	}
	
	public static void convert(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, File libFile, boolean writeBlib, boolean alignBetweenFiles) {
		ArrayList<File> featureFiles=new ArrayList<File>();
		SearchJobData representativeJob=null;
		for (int i=0; i<pecanJobs.size(); i++) {
			SearchJobData job=pecanJobs.get(i);
			if (!job.hasBeenRun()) {
				continue;
			}
			if (representativeJob==null) {
				representativeJob=job;
			}
			featureFiles.add(job.getFeatureFile());
		}
		
		if (representativeJob==null) return;
		Logger.logLine("Using "+representativeJob.getDiaFile().getName()+" to extract representative search parameters");

		String filename=libFile.getName();
		if (filename.lastIndexOf('.')>0) {
			filename=filename.substring(0, filename.lastIndexOf('.'));
		}
		File bigFeatureFile=new File(representativeJob.getFeatureFile().getParentFile(), filename+"_concatenated_features.txt");
		File bigPercolatorFile=new File(representativeJob.getFeatureFile().getParentFile(), filename+"_concatenated_results.txt");
		
		SearchParameters parameters=representativeJob.getParameters();
		float threshold=parameters.getEffectivePercolatorThreshold();
		try {
			ArrayList<PercolatorPeptide> passingPeptides;
			if (featureFiles.size()==1) {
				// if there's only one file then don't need to re-run percolator
				passingPeptides=PercolatorReader.getPassingPeptidesFromTSV(representativeJob.getOutputFile(), threshold);
			} else if (bigPercolatorFile.exists()&&bigPercolatorFile.canRead()) {
				// if we've already run percolator then don't need to re-run percolator
				passingPeptides=PercolatorReader.getPassingPeptidesFromTSV(bigPercolatorFile, threshold);
			} else {
				TableConcatenator.concatenateTables(featureFiles, bigFeatureFile);
				passingPeptides=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorLocation(), bigFeatureFile, bigPercolatorFile, threshold);
			}
			
			ArrayList<ProteinGroup> proteins=ParsimonyProteinGrouper.groupProteins(passingPeptides);
			Logger.logLine("Identified "+passingPeptides.size()+" peptides ("+proteins.size()+" proteins) across all files at a "+(threshold*100.0f)+"% FDR threshold.");

			Optional<PeakLocationInferrer> inferrer;
			if (alignBetweenFiles) {
			if (pecanJobs.size()>1) {
				Logger.logLine("Inferring peak boundaries across files...");
				inferrer=Optional.of(PeakLocationInferrer.getAlignmentData(new EmptyProgressIndicator(), pecanJobs, passingPeptides, parameters));
				Logger.logLine("...Finished peak inference.");
			} else {
				Logger.logLine("Only processing one file so no peak inference is necessary.");
				inferrer=Optional.empty();
			}
			} else {
				Logger.logLine("User requested no RT alignment between files.");
				inferrer=Optional.empty();
			}
			
			if (writeBlib) {
				convertBlib(progress, pecanJobs, libFile, Optional.of(passingPeptides), inferrer);
			} else {
				convertElib(progress, pecanJobs, libFile, Optional.of(passingPeptides), inferrer, proteins);
			}
			progress.update(passingPeptides.size()+" peptides identified at "+(threshold*100.0f)+"% FDR", 1.0f);
		} catch (IOException ioe) {
			Logger.errorLine("Error creating concatenated feature file");
			Logger.errorException(ioe);
		} catch (InterruptedException ie) {
			Logger.errorLine("Error creating concatenated feature file");
			Logger.errorException(ie);
		}
	}
	
	static void convertBlib(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, File blibFile, Optional<ArrayList<PercolatorPeptide>> passingPeptides, Optional<PeakLocationInferrer> inferrer) {
		try {
			BlibFile blib=new BlibFile();
			blib.openFile();
			blib.setUserFile(blibFile);
			blib.dropIndices();
			int[] counterTotals=new int[] {0,0,0};

			float increment=1.0f/pecanJobs.size();
			for (int i=0; i<pecanJobs.size(); i++) {
				SearchJobData job=pecanJobs.get(i);
				if (!job.hasBeenRun()) {
					continue;
				}
				ProgressIndicator subProgress=new SubProgressIndicator(progress, increment);
				
				ArrayList<PercolatorPeptide> globalPassingPeptides;
				ArrayList<PercolatorPeptide> localPassingPeptides=PercolatorReader.getPassingPeptidesFromTSV(job.getOutputFile(), pecanJobs.get(i).getParameters().getEffectivePercolatorThreshold());
				if (passingPeptides.isPresent()) {
					globalPassingPeptides=passingPeptides.get();
				} else {
					globalPassingPeptides=localPassingPeptides;
				}
				
				counterTotals=convertFileBlib(subProgress, job, globalPassingPeptides, localPassingPeptides, counterTotals, inferrer, blib);
			}

			blib.createIndices();
			blib.saveFile();
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
	 * @param subProgress
	 * @param job
	 * @param globalPassingPeptides
	 * @param localPassingPeptides
	 * @param counterTotals
	 * @param blib
	 * @param libraryFile
	 * @return
	 * @throws IOException
	 * @throws SQLException
	 */
	static int[] convertFileBlib(ProgressIndicator subProgress, SearchJobData job, ArrayList<PercolatorPeptide> globalPassingPeptides, ArrayList<PercolatorPeptide> localPassingPeptides, int[] counterTotals, Optional<PeakLocationInferrer> inferrer, BlibFile blib) throws IOException, SQLException {
		File diaFile=job.getDiaFile();
		Logger.logLine("Reading Percolator Results from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Reading Percolator Results", 0.0f);

		StripeFileInterface stripeFile=MzmlToDIAConverter.getFile(diaFile, job.getParameters());
		Logger.logLine("Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides", 0.00001f);

		LibraryInterface library=null;
		if (job instanceof EncyclopediaJobData) {
			library=((EncyclopediaJobData)job).getLibrary();
		}
		//ArrayList<IntegratedLibraryEntry> libraryEntries=SearchFeatureReader.parseSearchFeatures(featureFile, globalPassingPeptides, localPassingPeptides, stripeFile, Optional.ofNullable((LibraryFile)null), job.getParameters());
		ArrayList<IntegratedLibraryEntry> libraryEntries=PeptideQuantExtractor.parseSearchFeatures(subProgress, job, true, globalPassingPeptides, localPassingPeptides, inferrer, stripeFile, library, job.getParameters());
		stripeFile.close();
		
		File integrationFile=new File(diaFile.getAbsolutePath()+".integration.txt");

		PrintWriter writer=new PrintWriter(integrationFile, "UTF-8");
		writer.println("File Name\tPeptide Modified Sequence\tMin Start Time\tMax End Time\tPrecursor Charge\tPrecursorIsDecoy\tIon Count\tRetention Time Center\tTIC");
		
		for (IntegratedLibraryEntry entry : libraryEntries) {
			writer.println(diaFile.getName()+"\t"+entry.getPeptideModSeq()+"\t"+entry.getRtRange().getStart()/60f+"\t"+entry.getRtRange().getStop()/60f+"\t"+entry.getPrecursorCharge()+"\tFALSE\t"+entry.getIonCount()+"\t"+entry.getRetentionTime()/60f+"\t"+entry.getTIC());
		}
		writer.flush();
		writer.close();
		
		ArrayList<LibraryEntry> recasted=new ArrayList<LibraryEntry>();
		for (IntegratedLibraryEntry entry : libraryEntries) {
			recasted.add(entry);
		}

		Logger.logLine("Writing Skyline BLIB from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Writing Skyline BLIB", 0.99999f);

		counterTotals=blib.addLibrary(job, recasted, counterTotals[0], counterTotals[1], counterTotals[2]);
		subProgress.update(diaFile.getName()+": Finished writing to Skyline BLIB at"+new Date().toString(), 1.0f);
		return counterTotals;
	}
	
	static void convertElib(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, File elibFile, Optional<ArrayList<PercolatorPeptide>> passingPeptides, Optional<PeakLocationInferrer> inferrer, ArrayList<ProteinGroup> proteins) {
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
				ArrayList<PercolatorPeptide> localPassingPeptides=PercolatorReader.getPassingPeptidesFromTSV(job.getOutputFile(), pecanJobs.get(i).getParameters().getEffectivePercolatorThreshold());
				if (passingPeptides.isPresent()) {
					globalPassingPeptides=passingPeptides.get();
				} else {
					globalPassingPeptides=localPassingPeptides;
				}

				Logger.logLine(job.getDiaFile().getName()+": Number of global peptides: "+globalPassingPeptides.size()+" vs local peptides: "+localPassingPeptides.size());
				
				convertFileElib(subProgress, job, globalPassingPeptides, localPassingPeptides, inferrer, elib);
			}
			
			elib.setSources(pecanJobs);

			elib.createIndices();
			elib.saveAsFile(elibFile);
			
			if (pecanJobs.size()>1) {
				try {
					LibraryReportExtractor.extractMatrix(elib, proteins);
				} catch (DataFormatException e) {
					Logger.errorException(e);
				}
			}
			
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
	 * @param subProgress
	 * @param job
	 * @param globalPassingPeptides
	 * @param localPassingPeptides
	 * @param elib
	 * @param libraryFile
	 * @throws IOException
	 * @throws SQLException
	 */
	static void convertFileElib(ProgressIndicator subProgress, SearchJobData job, ArrayList<PercolatorPeptide> globalPassingPeptides, ArrayList<PercolatorPeptide> localPassingPeptides, Optional<PeakLocationInferrer> inferrer, LibraryFile elib) throws IOException, SQLException {
		File diaFile=job.getDiaFile();
		Logger.logLine("Reading Percolator Results from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Reading Percolator Results", 0.0f);

		StripeFileInterface stripeFile=MzmlToDIAConverter.getFile(diaFile, job.getParameters());
		Logger.logLine("Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides", 0.00001f);
		elib.addTIC(stripeFile);
		
		LibraryInterface library=null;
		if (job instanceof EncyclopediaJobData) {
			library=((EncyclopediaJobData)job).getLibrary();
		}
		ArrayList<IntegratedLibraryEntry> libraryEntries=PeptideQuantExtractor.parseSearchFeatures(subProgress, job, false, globalPassingPeptides, localPassingPeptides, inferrer, stripeFile, library, job.getParameters());
		stripeFile.close();
		
		Logger.logLine("Writing Encyclopedia ELIB from "+diaFile.getName()+" ("+libraryEntries.size()+" entries)...");
		subProgress.update(diaFile.getName()+": Writing Encyclopedia ELIB", 0.99999f);

		elib.addIntegratedEntries(libraryEntries, inferrer);
		Logger.logLine("Finished writing to Encyclopedia ELIB at"+new Date().toString());
		subProgress.update(diaFile.getName()+": Finished writing to Encyclopedia ELIB at"+new Date().toString(), 1.0f);
	}
}
