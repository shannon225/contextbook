package edu.washington.gs.maccoss.encyclopedia.cli;

import java.io.File;
import java.util.HashMap;
import java.util.StringTokenizer;

import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.ParsingUtils;
import edu.washington.gs.maccoss.encyclopedia.filewriters.LibraryUtilities;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class AdjustLibraryForPTMs {
	private static final boolean IS_FIXED=true;
	private static final boolean IS_COMBINED=false;
	
	public static void main(String[] args) {
		HashMap<String, String> arguments= CommandLineParser.parseArguments(args);
		if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("Adjust Library For PTMs ");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput library file");
			Logger.timelessLogLine("\t-o\toutput library file");
			Logger.timelessLogLine("\t-ptm\tcomma delimited modifications (for example \"K=8.014199,R=10.008269\")");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-fixed\t(default: "+IS_FIXED+")");
			Logger.timelessLogLine("\t-combineWithUnmodified\t(default: "+IS_COMBINED+")");
		} else {
			convert(arguments);
		}
	}

	public static void convert(HashMap<String, String> arguments) {
		if (!arguments.containsKey("-i")) {
			Logger.errorLine("You are required to specify an input library file (-i)");
			System.exit(1);
		}
		if (!arguments.containsKey("-o")) {
			Logger.errorLine("You are required to specify an output library file (-o)");
			System.exit(1);
		}
		if (!arguments.containsKey("-ptm")) {
			Logger.errorLine("You are required to specify modifications (-ptm), for example \"K=8.014199,R=10.008269\"");
			System.exit(1);
		}
		File inputFile=new File(arguments.get("-i"));
		File outputFile=new File(arguments.get("-o"));

		TCharDoubleHashMap ptms=new TCharDoubleHashMap();
		String value=arguments.get("-ptm");
		try {
			StringTokenizer st=new StringTokenizer(value, ",");
			while (st.hasMoreTokens()) {
				String token=st.nextToken();
				char aa=token.charAt(0);
				double mass=Double.parseDouble(token.substring(2)); // +1 for '=' (could actually be any deliminator)
				ptms.put(aa, mass);
			}
		} catch (Exception e) {
			throw new EncyclopediaException("Error parsing PTMs from ["+value+"]", e);
		}
		
		boolean isFixed=ParsingUtils.getBoolean("-fixed", arguments, IS_FIXED);
		boolean combineWithUnmodified=ParsingUtils.getBoolean("-combineWithUnmodified", arguments, IS_COMBINED);
			
		try {
			if (inputFile.exists()) {

				LibraryInterface library=BlibToLibraryConverter.getFile(inputFile);
				LibraryUtilities.modifyLibrary(outputFile, ptms, isFixed, combineWithUnmodified, library);
			} else {
				Logger.logLine("You must specify an existing ELIB or DLIB library file!");
			}
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}
}
