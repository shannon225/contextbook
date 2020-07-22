package edu.washington.gs.maccoss.encyclopedia.commandlineinterface;

import java.io.File;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;

import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaToPrositCSVParametersParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filewriters.StripeFileMerger;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class MergeDIAFiles {
	public static final String deliminator = ":";
	public static void main(String[] args) {
		HashMap<String, String> arguments= CommandLineParser.parseArguments(args);
		if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("CLI for Convert -> Merge DIA files");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .mzML files ("+deliminator+" deliminated)");
			Logger.timelessLogLine("\t-o\toutput .DIA file");

			TreeMap<String, String> defaults= new TreeMap<>(FastaToPrositCSVParametersParser.getDefaultParameters());
			int maxWidth=0;
			for (String key : defaults.keySet()) {
				if (key.length()>maxWidth) maxWidth=key.length();
			}
			maxWidth += 1;
			for (Entry<String, String> entry : defaults.entrySet()) {
				Logger.timelessLogLine("\t"+ General.formatCellToWidth(entry.getKey(), maxWidth)+" (default: "+entry.getValue()+")");
			}
		} else {
			convert(arguments);
		}
	}

	public static void convert(HashMap<String, String> arguments) {
		if (!arguments.containsKey("-i")) {
			Logger.errorLine("You are required to specify a FASTA file (-i)");
			System.exit(1);
		}

		File outputFile = new File(arguments.get("-o"));
		String[] inputs=arguments.get("-i").split(deliminator);
		File[] inputFiles=new File[inputs.length];
		
		for (int i = 0; i < inputFiles.length; i++) {
			inputFiles[i]=new File(inputs[i]);
			if (!inputFiles[i].exists()) {
				Logger.logLine("Input ["+inputs[i]+"] is not a file!");
				System.exit(1);
			}
		}

		try {
			StripeFileMerger.merge(inputFiles, outputFile, PecanParameterParser.parseParameters(PecanParameterParser.getDefaultParameters()));
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}
}
