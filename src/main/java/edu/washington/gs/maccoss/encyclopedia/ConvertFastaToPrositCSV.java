package edu.washington.gs.maccoss.encyclopedia;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaToPrositCSVParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaToPrositCSVParametersParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PrositCSVWriter;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

import java.io.File;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;

public class ConvertFastaToPrositCSV {
	public static void main(String[] args) {
		HashMap<String, String> arguments= CommandLineParser.parseArguments(args);
		if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("CLI for Convert -> Create Prosit CSV From FASTA");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .FASTA file");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-o\toutput .csv file");

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

		String outputFile = null;
		if (arguments.containsKey("-o")) {
			outputFile = arguments.get("-o");
		}

		File fastaFile = new File(arguments.get("-i"));
		FastaToPrositCSVParameters params= FastaToPrositCSVParametersParser.parseParameters(arguments);
		try {
			if (fastaFile.exists()) {
				PrositCSVWriter.writeCSV(outputFile, fastaFile, params.getEnzyme(), params.getDefaultNCE(), params.getDefaultCharge(),
						params.getMinCharge(), params.getMaxCharge(), params.getMaxMissedCleavage(),
						new Range(params.getMinMz(), params.getMaxMz()), false);
				Logger.logLine("Finished conversion from FASTA to Prosit CSV");
			} else {
				Logger.logLine("You must specify a FASTA file!");
			}
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}
}
