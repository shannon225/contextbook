package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;

public interface LibrarySearchJobData {
	public LibraryInterface getLibrary();
	public File getResultLibrary();
}
