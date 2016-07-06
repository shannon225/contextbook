package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;

public class EncyclopediaJobData extends SearchJobData {
	private final LibraryInterface library;
	private final LibraryScoringFactory taskFactory;
	
	public EncyclopediaJobData(File diaFile, LibraryInterface library, LibraryScoringFactory taskFactory) {
		super(diaFile, new File(diaFile.getAbsolutePath()+".features.txt"), new File(diaFile.getAbsolutePath()+".encyclopedia.txt"), taskFactory.getParameters(), taskFactory.getVersion());
		this.library=library;
		this.taskFactory=taskFactory;
	}

	public EncyclopediaJobData(File diaFile, LibraryInterface library, File featureFile, File outputFile, LibraryScoringFactory taskFactory) {
		super(diaFile, featureFile, outputFile, taskFactory.getParameters(), taskFactory.getVersion());
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
		return new File(getOutputFile().getAbsolutePath()+".first_round.txt");
	}
	
	public File getResultLibrary() {
		return new File(getOutputFile().getAbsolutePath()+".elib");
	}
}