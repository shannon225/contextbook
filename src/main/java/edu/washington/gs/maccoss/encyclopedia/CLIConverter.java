package edu.washington.gs.maccoss.encyclopedia;

import edu.washington.gs.maccoss.encyclopedia.cli.*;
import edu.washington.gs.maccoss.encyclopedia.utils.*;
import java.util.*;

public class CLIConverter {
	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		
		if (arguments.containsKey("-prositcsvtolibrary")||arguments.containsKey("-prositCSVToLibrary")) {
			ConvertPrositCSVToLibrary.main(args);
		} else if (arguments.containsKey("-libtoblib")||arguments.containsKey("-libraryToBlib")) {
			ConvertLibraryToBlib.main(args);
		} else if (arguments.containsKey("-fastatoprositcsv")||arguments.containsKey("-fastaToPrositCSV")) {
			ConvertFastaToPrositCSV.main(args);
		} else if (arguments.containsKey("-msptolib")||arguments.containsKey("-mspToLib")) {
			ConvertMSPToLibrary.main(args);
		} else if (arguments.containsKey("-blibtolib")||arguments.containsKey("-blibToLib")) {
			ConvertBLIBToLibrary.main(args);
		} else if (arguments.containsKey("-mergeDIA")||arguments.containsKey("-processDIA")) {
			PreprocessDIAFiles.main(args);
		} else if (arguments.containsKey("-adjustLibraryForPTMs")) {
			AdjustLibraryForPTMs.main(args);
		} else if (arguments.containsKey("-mergeLibraries")) {
			MergeLibraryFiles.main(args);
		} else if (arguments.containsKey("-openswathTSVToLibrary")) {
			ConvertOpenSwathToLibrary.main(args);
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("EncyclopeDIA Converter Help");
			Logger.timelessLogLine("EncyclopeDIA Converter allows to convert files from one format to another.");
			Logger.timelessLogLine("Available converters: ");
			Logger.timelessLogLine("\t-prositCSVToLibrary\tConvert Prosit/Spectronaut CSV to library (use -convert -prositcsvtolibrary -h for help)");
			Logger.timelessLogLine("\t-blibToLib\tConvert Convert BLIB to Library (use -convert -blibToLib -h for help)");
			Logger.timelessLogLine("\t-mspToLib\tConvert Convert SPTXT/MSP to Library (use -convert -mspToLib -h for help)");
			Logger.timelessLogLine("\t-openswathTSVToLibrary\tConvert OpenSwath TSV to library (use -convert -openswathTSVToLibrary -h for help)");
			Logger.timelessLogLine("\t-libraryToBlib\tConvert library to BLIB (use -convert -libtoblib -h for help)");
			Logger.timelessLogLine("\t-mergeLibraries\tMerge multiple DLIB libraries into a single DLIB (use -convert -mergeLibraries -h for help)");
			Logger.timelessLogLine("\t-fastaToPrositCSV\tConvert FASTA to Prosit CSV (use -convert -fastatoprositcsv -h for help)");
			Logger.timelessLogLine("\t-processDIA\tPreprocess .MZMLs or merge .MZML or .DIA gas-phase fractions (use -convert -processDIA -h for help)");
			Logger.timelessLogLine("\t-adjustLibraryForPTMs\tadd PTMs (or SILAC masses) to library (use -convert -adjustLibraryForPTMs -h for help)");
			System.exit(1);
		}
	}
}
