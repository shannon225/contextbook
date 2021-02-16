package edu.washington.gs.maccoss.encyclopedia;

import edu.washington.gs.maccoss.encyclopedia.cli.AdjustLibraryForPTMs;
import edu.washington.gs.maccoss.encyclopedia.cli.ConvertFastaToPrositCSV;
import edu.washington.gs.maccoss.encyclopedia.cli.ConvertLibraryToBlib;
import edu.washington.gs.maccoss.encyclopedia.cli.ConvertOpenSwathToLibrary;
import edu.washington.gs.maccoss.encyclopedia.cli.ConvertPrositCSVToLibrary;
import edu.washington.gs.maccoss.encyclopedia.cli.PreprocessDIAFiles;
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
		} else if (arguments.containsKey("-mergeDIA")||arguments.containsKey("-processDIA")) {
			PreprocessDIAFiles.main(args);
		} else if (arguments.containsKey("-adjustLibraryForPTMs")) {
			AdjustLibraryForPTMs.main(args);
		} else if (arguments.containsKey("-openswathTSVToLibrary")) {
			ConvertOpenSwathToLibrary.main(args);
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("EncyclopeDIA Converter Help");
			Logger.timelessLogLine("EncyclopeDIA Converter allows to convert files from one format to another.");
			Logger.timelessLogLine("Available converters: ");
			Logger.timelessLogLine("\t-prositCSVToLibrary\tConvert Prosit/Spectronaut CSV to library (use -convert -prositcsvtolibrary -h for help)");
			Logger.timelessLogLine("\t-openswathTSVToLibrary\tConvert OpenSwath TSV to library (use -convert -openswathTSVToLibrary -h for help)");
			Logger.timelessLogLine("\t-libraryToBlib\tConvert library to BLIB (use -convert -libtoblib -h for help)");
			Logger.timelessLogLine("\t-fastaToPrositCSV\tConvert FASTA to Prosit CSV (use -convert -fastatoprositcsv -h for help)");
			Logger.timelessLogLine("\t-processDIA\tPreprocess .MZMLs or merge .MZML or .DIA gas-phase fractions (use -convert -processDIA -h for help)");
			Logger.timelessLogLine("\t-adjustLibraryForPTMs\tadd PTMs (or SILAC masses) to library (use -convert -adjustLibraryForPTMs -h for help)");
			System.exit(1);
		}
	}
}
