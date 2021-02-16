package edu.washington.gs.maccoss.encyclopedia.cli;

import java.io.File;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.OpenSwathTSVToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class ConvertOpenSwathToLibrary {
	public static void main(String[] args) {
		HashMap<String, String> arguments= CommandLineParser.parseArguments(args);
		if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("Convert OpenSwath TSV to library");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .TSV file");
			Logger.timelessLogLine("\t-f\ttaxon .FASTA database");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-o\toutput .dlib file");
		} else {
			convert(arguments);
		}
	}

	public static void convert(HashMap<String, String> arguments) {
		if (!arguments.containsKey("-i") || !arguments.containsKey("-f")) {
			Logger.errorLine("You are required to specify an input TSV file (-i) and a fasta file (-f)");
			System.exit(1);
		}

		File inputFile = new File(arguments.get("-i"));
		File fastaFile = new File(arguments.get("-f"));
		File outputFile;
		if (arguments.containsKey("-o")) {
			outputFile = new File(arguments.get("-o"));
		} else {
			String absolutePath=inputFile.getAbsolutePath();
			outputFile = new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		}

		SearchParameters parameters= SearchParameterParser.parseParameters(arguments);
		try {
			OpenSwathTSVToLibraryConverter.convertFromOpenSwathTSV(inputFile, fastaFile, outputFile, parameters);
			Logger.logLine("Finished reading " + inputFile.getName());
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}
}
