package edu.washington.gs.maccoss.encyclopedia.cli;

import java.io.File;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filewriters.StripeFileMerger;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class PreprocessDIAFiles {
	public static final String deliminator = ":";
	public static void main(String[] args) {
		HashMap<String, String> arguments= CommandLineParser.parseArguments(args);
		if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("CLI for Convert -> Preprocess DIA files");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .mzML or .DIA files (merge multiple files by deliminating with \""+deliminator+"\")");
			Logger.timelessLogLine("\t-o\toutput merged .DIA file (only used for merging multiple files)");
		} else {
			convert(arguments);
		}
	}

	public static void convert(HashMap<String, String> arguments) {
		if (!arguments.containsKey("-i")) {
			Logger.errorLine("You are required to specify a FASTA file (-i)");
			System.exit(1);
		}

		String[] inputs=arguments.get("-i").split(deliminator);
		File[] inputFiles=new File[inputs.length];
		
		for (int i = 0; i < inputFiles.length; i++) {
			inputFiles[i]=new File(inputs[i]);
			if (!inputFiles[i].exists()) {
				Logger.logLine("Input ["+inputs[i]+"] is not a file!");
				System.exit(1);
			}
		}

		File outputFile = null;
		if (inputFiles.length==1) {
			// expect no output file (will write in place)
			if (arguments.get("-o")!=null) {
				Logger.logLine("When using one input, do not specify an output file!");
				System.exit(1);
			}
		} else {
			outputFile=new File(arguments.get("-o"));
		}
		
		SearchParameters parameters=SearchParameterParser.parseParameters(arguments);

		try {
			if (inputFiles.length==1) {
				StripeFileGenerator.getFile(inputFiles[0], parameters);
			} else {
				StripeFileMerger.merge(inputFiles, outputFile, parameters);
			}
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}
}
