package edu.washington.gs.maccoss.encyclopedia.cli;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaToPrositCSVParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaToPrositCSVParametersParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaLibraryPredictionClient;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;

public class ConvertFastaToKoinaPrositLibrary {
	public static void main(String[] args) {
		HashMap<String, String> arguments= CommandLineParser.parseArguments(args);
		if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("CLI for Convert -> Create Prosit Library From FASTA");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .FASTA file");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-o\toutput .dlib file");

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
		ArrayList<KoinaFeaturePredictionModel> models = KoinaLibraryPredictionClient.getDefaultModels();

		String outputFile = null;
		if (arguments.containsKey("-o")) {
			outputFile = arguments.get("-o");
		}

		File fastaFile = new File(arguments.get("-i"));
		FastaToPrositCSVParameters params= FastaToPrositCSVParametersParser.parseParameters(arguments);
		try {
			if (fastaFile.exists()) {
				KoinaLibraryPredictionClient.writeLibrary(models, outputFile, fastaFile, params.getEnzyme(), params.getDefaultNCE(), params.getDefaultCharge(),
						params.getMinCharge(), params.getMaxCharge(), params.getMaxMissedCleavage(), new Range(params.getMinMz(), params.getMaxMz()), 
						params.isAdjustNCEForDIA(), params.isAddDecoys(), SearchParameterParser.getDefaultParametersObject(), new EmptyProgressIndicator(false));
				
				Logger.logLine("Finished conversion from FASTA to Prosit dlib");
			} else {
				Logger.logLine("You must specify a FASTA file!");
			}
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}
}
