package edu.washington.gs.maccoss.encyclopedia.cli;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.ParsingUtils;
import edu.washington.gs.maccoss.encyclopedia.filewriters.LibraryUtilities;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;

public class MergeLibraryFiles {
	public static final String deliminator = ":";
	private static final boolean IS_RT_ALIGN=false;
	private static final boolean IS_HIGHER_SCORES_ARE_BETTER=false;
	private static final boolean IS_REMOVE_DUPLICATES=true;
	
	public static void main(String[] args) {
		HashMap<String, String> arguments= CommandLineParser.parseArguments(args);
		if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("Merge Library Files");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput library files or directories (merge multiple files by deliminating with \""+deliminator+"\")");
			Logger.timelessLogLine("Other Parameters: ");
			Logger.timelessLogLine("\t-o\toutput .DLIB file");
			Logger.timelessLogLine("\t-rtAlign\t(default: "+IS_RT_ALIGN+")");
			Logger.timelessLogLine("\t-removeDuplicates\t(default: "+IS_REMOVE_DUPLICATES+")");
			Logger.timelessLogLine("\t-higherScoresAreBetter\t(default: "+IS_HIGHER_SCORES_ARE_BETTER+"), only used if removeDuplicates=true");
		} else {
			convert(arguments);
		}
	}

	public static void convert(HashMap<String, String> arguments) {
		if (!arguments.containsKey("-i")) {
			Logger.errorLine("You are required to specify an input library file (-i)");
			System.exit(1);
		}
		
		String[] inputs=arguments.get("-i").split(deliminator);
		ArrayList<File> files=new ArrayList<>();
		File firstFile=null;
		
		SimpleFilenameFilter libFilter=new SimpleFilenameFilter(LibraryFile.DLIB, LibraryFile.ELIB);
		for (int i = 0; i < inputs.length; i++) {
			File f=new File(inputs[i]);
			if (firstFile==null) {
				firstFile=f;
			}
			if (f.isDirectory()) {
				File[] fs=f.listFiles();
				for (File file : fs) {
					if (libFilter.accept(file.getName())) {
						files.add(file);
					}
				}
			}
			if (libFilter.accept(f.getName())) {
				files.add(f);
			}
		}

		File saveFile;
		if (arguments.containsKey("-o")) {
			saveFile = new File(arguments.get("-o"));
		} else {
			String absolutePath;
			if (firstFile.isDirectory()) {
				absolutePath=firstFile.getParent()+File.pathSeparator+firstFile.getName();
			} else {
				String path=firstFile.getAbsolutePath();
				absolutePath=path.substring(0, path.lastIndexOf('.'));
			}
			saveFile=new File(absolutePath+"_merged"+LibraryFile.DLIB);
		}
		

		boolean rtAlign=ParsingUtils.getBoolean("-rtAlign", arguments, IS_RT_ALIGN);

		boolean removeDuplicates=ParsingUtils.getBoolean("-removeDuplicates", arguments, IS_REMOVE_DUPLICATES);
		boolean higherScoresAreBetter=ParsingUtils.getBoolean("-higherScoresAreBetter", arguments, IS_HIGHER_SCORES_ARE_BETTER);

		try {
			Logger.logLine("Merging "+files.size()+" libraries into "+saveFile.getName());
			LibraryUtilities.mergeLibraries(new EmptyProgressIndicator(false), files, saveFile, rtAlign, removeDuplicates, higherScoresAreBetter);
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}
}
