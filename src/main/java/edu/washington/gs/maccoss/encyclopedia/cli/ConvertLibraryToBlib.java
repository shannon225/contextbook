package edu.washington.gs.maccoss.encyclopedia.cli;

import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryToBlibConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

import java.io.File;
import java.util.HashMap;

public class ConvertLibraryToBlib {
	public static void main(String[] args) {
		HashMap<String, String> arguments= CommandLineParser.parseArguments(args);
		if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("Convert Library to BLIB ");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput library file");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-o\toutput .blib file");
		} else {
			convert(arguments);
		}
	}

	public static void convert(HashMap<String, String> arguments) {
		if (!arguments.containsKey("-i")) {
			Logger.errorLine("You are required to specify an input library file (-i)");
			System.exit(1);
		}
		File elibFile=new File(arguments.get("-i"));

		File blibFile;
		if (arguments.containsKey("-o")) {
			blibFile = new File(arguments.get("-o"));
		} else {
			String absolutePath=elibFile.getAbsolutePath();
			blibFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+ BlibFile.BLIB);
		}

		try {
			if (elibFile.exists()) {
				LibraryToBlibConverter.convert(elibFile, blibFile);
			} else {
				Logger.logLine("You must specify an existing ELIB or DLIB library file!");
			}
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}
}
