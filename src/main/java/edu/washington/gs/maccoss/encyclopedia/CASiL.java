package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.ExecutionException;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoEncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.FileLogRecorder;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class CASiL {

	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		if (arguments.size()==0) {
			SearchGUIMain.runGUI(ProgramType.CASiL);
			
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("CASiL Help");
			Logger.timelessLogLine("CASiL is a Chromatogram Aligned Site Localizing search engine for DIA data.");
			Logger.timelessLogLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .DIA or .MZML file");
			Logger.timelessLogLine("\t-l\tlibrary .ELIB file");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-o\toutput report file (default: [input file]"+EncyclopediaJobData.OUTPUT_FILE_SUFFIX+")");
			
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
			Logger.logLine("CASiL version "+PhosphoEncyclopediaOneScoringFactory.version);
			System.exit(1);
			
		} else {
			if (!arguments.containsKey(Encyclopedia.INPUT_DIA_TAG)||!arguments.containsKey(Encyclopedia.TARGET_LIBRARY_TAG)) {
				Logger.errorLine("You are required to specify an input file ("+Encyclopedia.INPUT_DIA_TAG+") and a library file ("+Encyclopedia.TARGET_LIBRARY_TAG+")");
				System.exit(1);
			}

			File diaFile=new File(arguments.get(Encyclopedia.INPUT_DIA_TAG));
			File libraryFile=new File(arguments.get(Encyclopedia.TARGET_LIBRARY_TAG));

			File outputFile;
			if (arguments.containsKey(Encyclopedia.OUTPUT_RESULT_TAG)) {
				outputFile=new File(arguments.get(Encyclopedia.OUTPUT_RESULT_TAG));
			} else {
				outputFile=new File(diaFile.getAbsolutePath()+EncyclopediaJobData.OUTPUT_FILE_SUFFIX);
			}

			try {
				FileLogRecorder logRecorder=new FileLogRecorder(new File(outputFile.getAbsolutePath()+EncyclopediaJobData.LOG_FILE_SUFFIX));
				Logger.addRecorder(logRecorder);
	
				SearchParameters parameters=SearchParameterParser.parseParameters(arguments);
				if (!parameters.getLocalizingModification().isPresent()) {
					Logger.errorLine("You are required to specify one localization modification ("+PeptideModification.getShortnameList()+")");
					System.exit(1);
				}

				Logger.logLine("Setting up localization engne...");
				StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
				PhosphoLocalizer localizer=new PhosphoLocalizer(stripefile, parameters.getLocalizingModification().get(), parameters);
				LibraryScoringFactory factory=new PhosphoEncyclopediaOneScoringFactory(parameters, localizer);
				
				Logger.logLine("CASiL version "+factory.getVersion());
	
				Logger.logLine("Parameters:");
				Logger.logLine(" "+Encyclopedia.INPUT_DIA_TAG+" "+diaFile.getAbsolutePath());
				Logger.logLine(" "+Encyclopedia.TARGET_LIBRARY_TAG+" "+libraryFile.getAbsolutePath());
				Logger.logLine(" "+Encyclopedia.OUTPUT_RESULT_TAG+" "+outputFile.getAbsolutePath());
				Logger.logLine(parameters.toString());

				LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
				EncyclopediaJobData job=new EncyclopediaJobData(diaFile, library, outputFile, factory);
				runSearch(new EmptyProgressIndicator(), job);
			} catch (Exception e) {
				Logger.errorLine("Encountered Fatal Error!");
				Logger.errorException(e);
			} finally {
				Logger.close();
			}
		}
	}

	public static void runSearch(ProgressIndicator progress, EncyclopediaJobData job) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		if (!(job.getTaskFactory() instanceof PhosphoEncyclopediaOneScoringFactory)) {
			if (!job.getParameters().getLocalizingModification().isPresent()) {
				throw new EncyclopediaException("You are required to specify one localization modification ("+PeptideModification.getShortnameList()+")");
			}
			
			Logger.logLine("Setting up localization engne...");
			StripeFileInterface stripefile=StripeFileGenerator.getFile(job.getDiaFile(), job.getParameters());
			PhosphoLocalizer localizer=new PhosphoLocalizer(stripefile, job.getParameters().getLocalizingModification().get(), job.getParameters());
			LibraryScoringFactory factory=new PhosphoEncyclopediaOneScoringFactory(job.getParameters(), localizer);
			job=job.updateTaskFactory(factory);
		}
		
		File outputFile=job.getOutputFile();
		if (outputFile.exists()&&outputFile.canRead()) {
			try {
				ArrayList<PercolatorPeptide> passingPeptidesFromTSV=PercolatorReader.getPassingPeptidesFromTSV(outputFile, job.getParameters().getEffectivePercolatorThreshold());
				
				File elibFile=job.getResultLibrary();
				if (!elibFile.exists()) {
					progress.update("Writing elib result library...");
					Logger.logLine("Writing elib result library...");
					ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
					jobs.add(job);
					SearchToBLIB.convert(progress, jobs, elibFile, false, true);
				}
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
		
		SearchParameters parameters=job.getParameters();
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		if (parameters.isDDA()) {
			EncyclopediaDDA.runSearch(progress, job, stripefile);
		} else {
			Encyclopedia.runSearch(progress, job, stripefile);
		}
		stripefile.close();
	}
}
