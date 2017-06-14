package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoEncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.FileLogRecorder;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;

public class CAPSiL {

	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		if (arguments.size()==0) {
			SearchGUIMain.runGUI(ProgramType.CAPSiL);
			
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("CAPSiL Help");
			Logger.timelessLogLine("CAPSiL is a Chromatogram Aligned Phospho-Site Localizing search engine for DIA data.");
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
			Logger.logLine("CAPSiL version "+PhosphoEncyclopediaOneScoringFactory.version);
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

			arguments.put("-runPhosphoLocalization", "true");

			try {
				FileLogRecorder logRecorder=new FileLogRecorder(new File(outputFile.getAbsolutePath()+EncyclopediaJobData.LOG_FILE_SUFFIX));
				Logger.addRecorder(logRecorder);
	
				SearchParameters parameters=SearchParameterParser.parseParameters(arguments);

				Logger.logLine("Setting up localization engne...");
				StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
				PhosphoLocalizer localizer=new PhosphoLocalizer(stripefile, parameters);
				LibraryScoringFactory factory=new PhosphoEncyclopediaOneScoringFactory(parameters, localizer);
				
				Logger.logLine("CAPSiL version "+factory.getVersion());
	
				Logger.logLine("Parameters:");
				Logger.logLine(" "+Encyclopedia.INPUT_DIA_TAG+" "+diaFile.getAbsolutePath());
				Logger.logLine(" "+Encyclopedia.TARGET_LIBRARY_TAG+" "+libraryFile.getAbsolutePath());
				Logger.logLine(" "+Encyclopedia.OUTPUT_RESULT_TAG+" "+outputFile.getAbsolutePath());
				Logger.logLine(parameters.toString());

				LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
				EncyclopediaJobData job=new EncyclopediaJobData(diaFile, library, outputFile, factory);
				Encyclopedia.runSearch(new EmptyProgressIndicator(), job);
			} catch (Exception e) {
				Logger.errorLine("Encountered Fatal Error!");
				Logger.errorException(e);
			} finally {
				Logger.close();
			}
		}
	}
}
