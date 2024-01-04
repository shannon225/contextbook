package edu.washington.gs.maccoss.encyclopedia.cli;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ImmutableMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.OpenSwathTSVToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.StringUtils;

public class ConvertOpenSwathToLibrary {
	public static void main(String[] args) {
		Logger.logLine("Convert OpenSwath TSV to library");

		HashMap<String, String> arguments= CommandLineParser.parseArguments(args);
		if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .TSV file");
			Logger.timelessLogLine("\t-f\ttaxon .FASTA database");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-o\toutput .dlib file");
			Logger.timelessLogLine("");
			Logger.timelessLogLine("Input TSV format:");

			StringBuilder b = new StringBuilder();
			for (Map.Entry<String, List<String>> e :
					ImmutableMap.<String, List<String>>builder()
					.put("Peptide sequence and modifications", OpenSwathTSVToLibraryConverter.PEPTIDE_HEADERS)
					.put("Precursor charge", OpenSwathTSVToLibraryConverter.CHARGE_HEADERS)
					.put("Fragment m/z", OpenSwathTSVToLibraryConverter.FRAG_MZ_HEADERS)
					.put("Fragment intensity", OpenSwathTSVToLibraryConverter.FRAG_INTEN_HEADERS)
					.put("Retention time", OpenSwathTSVToLibraryConverter.RT_HEADERS)
					.put("Unique integer ID for each entry (optional)", OpenSwathTSVToLibraryConverter.TRANSITION_GROUP_HEADERS)
					.build().entrySet()
			) {
				String desc = e.getKey();
				List<String> aliases = e.getValue();

				if (aliases.isEmpty()) { // Should never happen
					continue;
				}

				b.setLength(0);
				b.append("\t");
				b.append(aliases.get(0));
				b.append(": ");
				b.append(desc);
				if (aliases.size() > 1) {
					b.append(" (aliases: ");
					for (int i = 1; i < aliases.size(); i++) {
						if (i > 1) {
							b.append(", ");
						}
						b.append(aliases.get(i));
					}
					b.append(")");
				}

				Logger.timelessLogLine(b.toString());
			}
		} else {
			convert(arguments);
		}
	}

	public static void convert(HashMap<String, String> arguments) {
		if (!arguments.containsKey("-i") || !arguments.containsKey("-f")) {
			Logger.errorLine("You are required to specify an input TSV file (-i) and a fasta file (-f)");
			System.exit(1);
		}

		File inputFile = new File(arguments.get("-i"));
		File fastaFile = new File(arguments.get("-f"));
		File outputFile;
		if (arguments.containsKey("-o")) {
			outputFile = new File(arguments.get("-o"));
		} else {
			String absolutePath=inputFile.getAbsolutePath();
			outputFile = new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		}

		SearchParameters parameters= SearchParameterParser.parseParameters(arguments);
		try {
			OpenSwathTSVToLibraryConverter.convertFromOpenSwathTSV(inputFile, fastaFile, outputFile, parameters);
			Logger.logLine("Finished reading " + inputFile.getName());
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}
}
