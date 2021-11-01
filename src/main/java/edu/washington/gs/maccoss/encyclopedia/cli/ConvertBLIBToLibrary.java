package edu.washington.gs.maccoss.encyclopedia.cli;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MSPReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

import java.io.File;
import java.util.HashMap;
import java.util.Optional;

public class ConvertBLIBToLibrary {
	public static void main(String[] args) {
		HashMap<String, String> arguments= CommandLineParser.parseArguments(args);
		if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("Convert Convert BLIB to Library");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .BLIB file");
			Logger.timelessLogLine("\t-f\ttaxon .FASTA database");
			Logger.timelessLogLine("\t-higherScoresAreBetter\t(true/false, default: [true])");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-o\toutput .dlib file");
			Logger.timelessLogLine("\t-r\tinput .IRTDB file");
		} else {
			convert(arguments);
		}
	}

	public static void convert(HashMap<String, String> arguments) {
		if (!arguments.containsKey("-i") || !arguments.containsKey("-f")) {
			Logger.errorLine("You are required to specify an input MSP file (-i) and a fasta file (-f)");
			System.exit(1);
		}

		File blibFile = new File(arguments.get("-i"));
		File fastaFile = new File(arguments.get("-f"));
		String higherScoresAreBetterString=arguments.get("-higherScoresAreBetter");
		boolean higherScoresAreBetter=true;
		if (higherScoresAreBetterString!=null) {
			if ("false".equalsIgnoreCase(higherScoresAreBetterString)||"f".equalsIgnoreCase(higherScoresAreBetterString)) {
				higherScoresAreBetter=false;
			}
		}
		
		File outputFile;
		if (arguments.containsKey("-o")) {
			outputFile = new File(arguments.get("-o"));
		} else {
			String absolutePath=blibFile.getAbsolutePath();
			outputFile = new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		}
		Optional<File> irtdbFile;
		if (arguments.containsKey("-r")) {
			irtdbFile=Optional.of(new File(arguments.get("-r")));
		} else {
			irtdbFile=Optional.empty();
		}

		SearchParameters parameters= SearchParameterParser.parseParameters(arguments);
		try {
			BlibToLibraryConverter.convert(blibFile, outputFile, irtdbFile, fastaFile, higherScoresAreBetter, parameters);
			Logger.logLine("Finished reading " + blibFile.getName());
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}
}
