package edu.washington.gs.maccoss.encyclopedia;

import edu.washington.gs.maccoss.encyclopedia.utils.*;
import java.util.*;

public class CLIConverter {
	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		
		if (arguments.containsKey("-prositcsvtolibrary")) {
			ConvertPrositCSVToLibrary.main(args);
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("EncyclopeDIA Converter Help");
			Logger.timelessLogLine("EncyclopeDIA Converter allows to convert files from one format to another.");
			Logger.timelessLogLine("Available converters: ");
			Logger.timelessLogLine("\t-prositcsvtolibrary\tConvert Prosit/Spectronaut CSV to library (use -convert -prositcsvtolibrary -h for help)");
			System.exit(1);
		}
	}
}
