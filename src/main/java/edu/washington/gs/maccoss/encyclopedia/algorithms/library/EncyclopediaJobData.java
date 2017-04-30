package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;

public class EncyclopediaJobData extends SearchJobData {
	public static final String LOG_FILE_SUFFIX=".log";
	public static final String OUTPUT_FILE_SUFFIX=".encyclopedia.txt";
	public static final String FEATURE_FILE_SUFFIX=".features.txt";
	private final LibraryInterface library;
	private final LibraryScoringFactory taskFactory;
	
	public EncyclopediaJobData(File diaFile, LibraryInterface library, LibraryScoringFactory taskFactory) {
		super(diaFile, new File(diaFile.getAbsolutePath()+FEATURE_FILE_SUFFIX), new File(diaFile.getAbsolutePath()+OUTPUT_FILE_SUFFIX), taskFactory.getParameters(), taskFactory.getVersion());
		this.library=library;
		this.taskFactory=taskFactory;
	}

	public EncyclopediaJobData(File diaFile, LibraryInterface library, File outputFile, LibraryScoringFactory taskFactory) {
		super(diaFile, new File(getOutputAbsolutePathPrefix(outputFile.getAbsolutePath())+FEATURE_FILE_SUFFIX), outputFile, taskFactory.getParameters(), taskFactory.getVersion());
		this.library=library;
		this.taskFactory=taskFactory;
	}
	
	public LibraryInterface getLibrary() {
		return library;
	}

	public LibraryScoringFactory getTaskFactory() {
		return taskFactory;
	}
	
	public File getFirstPassPercolator() {
		String absolutePath=getOutputAbsolutePathPrefix(getOutputFile().getAbsolutePath());
		return new File(absolutePath+".first_round.txt");
	}
	
	public File getResultLibrary() {
		String absolutePath=getOutputAbsolutePathPrefix(getOutputFile().getAbsolutePath());
		return new File(absolutePath+LibraryFile.ELIB);
	}

	public static String getOutputAbsolutePathPrefix(String absolutePath) {
		if (absolutePath.endsWith(OUTPUT_FILE_SUFFIX)) {
			absolutePath=absolutePath.substring(0, absolutePath.length()-OUTPUT_FILE_SUFFIX.length());
		}
		return absolutePath;
	}
	
	@Override
	public String getSearchType() {
		return "EncyclopeDIA";
	}
}