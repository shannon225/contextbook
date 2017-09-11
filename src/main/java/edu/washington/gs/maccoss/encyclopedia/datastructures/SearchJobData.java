package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;

/**
 * @author Seth Just
 * @since 2:50 PM 9/11/17
 */
public interface SearchJobData {
	File getDiaFile();

	File getFeatureFile();

	File getOutputFile();

	File getOutputDecoyFile();

	SearchParameters getParameters();

	String getVersion();

	boolean hasBeenRun();

	String getSearchType();
}
