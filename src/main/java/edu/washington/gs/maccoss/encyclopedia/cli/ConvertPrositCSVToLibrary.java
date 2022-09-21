package edu.washington.gs.maccoss.encyclopedia.cli;

import java.io.File;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SpectronautCSVToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class ConvertPrositCSVToLibrary {
	public static void main(String[] args) {
		HashMap<String, String> arguments= CommandLineParser.parseArguments(args);
		if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("Convert Prosit/Spectronaut/DIA-NN CSV to library");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .CSV file");
			Logger.timelessLogLine("\t-f\ttaxon .FASTA database");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-o\toutput .dlib file");
		} else {
			convert(arguments);
		}
	}

	public static void convert(HashMap<String, String> arguments) {
		if (!arguments.containsKey("-i") || !arguments.containsKey("-f")) {
			Logger.errorLine("You are required to specify an input CSV file (-i) and a fasta file (-f)");
			System.exit(1);
		}

		File csvFile = new File(arguments.get("-i"));
		File fastaFile = new File(arguments.get("-f"));
		File outputFile;
		if (arguments.containsKey("-o")) {
			outputFile = new File(arguments.get("-o"));
		} else {
			String absolutePath=csvFile.getAbsolutePath();
			outputFile = new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		}

		SearchParameters parameters= SearchParameterParser.parseParameters(arguments);
		try {
			SpectronautCSVToLibraryConverter.convertFromSpectronautCSV(csvFile, fastaFile, outputFile, parameters);
			Logger.logLine("Finished reading " + csvFile.getName());
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}
}
