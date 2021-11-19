package edu.washington.gs.maccoss.encyclopedia;

import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class Main {

	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		if (arguments.size()==0) {
			SearchGUIMain.runGUI(ProgramType.Global);
			
		} else if (arguments.size()==1&&arguments.containsKey(SearchParameters.ENABLE_ADVANCED_OPTIONS)) {
			SearchGUIMain.runGUI(ProgramType.Global, true);
			
		} else if (arguments.containsKey("-encyclopedia")) {
			Encyclopedia.main(args);
			
		} else if (arguments.containsKey("-thesaurus")) {
			Thesaurus.main(args);
			
		//} else if (arguments.containsKey("-scribe")) {
		//	Scribe.main(args);
			
		} else if (arguments.containsKey("-browser")) {
			DIABrowser.main(args);
		
		} else if (arguments.containsKey("-libexport")) {
			SearchToBLIB.main(args);

		} else if (arguments.containsKey("-convert")) {
			CLIConverter.main(args);
			
		} else if (arguments.containsKey("-pecan")) {
			Pecanpie.main(args);
			
		} else if (arguments.containsKey("-walnut")) {
			Pecanpie.main(args);
			
		} else if (arguments.containsKey("-xcordia")) {
			XCorDIA.main(args);
			
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("EncyclopeDIA Help");
			Logger.timelessLogLine("EncyclopeDIA is a library search engine for DIA data.");
			Logger.timelessLogLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .DIA or .MZML file");
			Logger.timelessLogLine("\t-l\tlibrary .ELIB file");
			Logger.timelessLogLine("Other Programs: ");
			Logger.timelessLogLine("\t-walnut\trun Walnut (use -walnut -h for Walnut help)");
			Logger.timelessLogLine("\t-thesaurus\trun Thesaurus (use -thesaurus -h for Thesaurus help)");
			Logger.timelessLogLine("\t-xcordia\trun XCorDIA (use -xcordia -h for XCorDIA help)");
			//Logger.timelessLogLine("\t-scribe\trun Scribe (use -scribe -h for Scribe help)");
			Logger.timelessLogLine("\t-browser\trun ELIB Browser (use -browser -h for ELIB Browser help)");
			Logger.timelessLogLine("\t-libexport\trun Library Export (use -libexport -h for Library Export help)");
			Logger.timelessLogLine("\t-convert\trun files converter (use -convert -h for help)");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-o\toutput report file (default: [input file]"+EncyclopediaJobData.OUTPUT_FILE_SUFFIX+")");

			TreeMap<String, String> defaults=new TreeMap<String, String>(SearchParameterParser.getDefaultParameters());
			int maxWidth=0;
			for (String key : defaults.keySet()) {
				if (key.length()>maxWidth) maxWidth=key.length();
			}
			for (Entry<String, String> entry : defaults.entrySet()) {
				Logger.timelessLogLine("\t"+General.formatCellToWidth(entry.getKey(), maxWidth)+" (default: "+entry.getValue()+")");
			}
			System.exit(1);
			
		} else {
			Encyclopedia.main(args);
		}
	}
}
