package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;

public class EncyclopediaJobData extends SearchJobData {
	private final File libraryFile;
	private final LibraryScoringFactory taskFactory;
	
	public EncyclopediaJobData(File diaFile, File libraryFile, LibraryScoringFactory taskFactory) {
		super(diaFile, new File(diaFile.getAbsolutePath()+".features.txt"), new File(diaFile.getAbsolutePath()+".percolator.txt"), taskFactory.getParameters(), taskFactory.getVersion());
		this.libraryFile=libraryFile;
		this.taskFactory=taskFactory;
	}

	public EncyclopediaJobData(File diaFile, File libraryFile, File featureFile, File outputFile, LibraryScoringFactory taskFactory) {
		super(diaFile, featureFile, outputFile, taskFactory.getParameters(), taskFactory.getVersion());
		this.libraryFile=libraryFile;
		this.taskFactory=taskFactory;
	}
	
	public File getLibraryFile() {
		return libraryFile;
	}

	public LibraryScoringFactory getTaskFactory() {
		return taskFactory;
	}
}