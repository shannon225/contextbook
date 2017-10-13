package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;

import java.io.File;

/**
 * @author Seth Just
 * @since 2:50 PM 9/11/17
 */
public interface SearchJobData {
	File getDiaFile();

	StripeFileInterface getDiaFileReader();

	PercolatorExecutionData getPercolatorFiles();

	SearchParameters getParameters();

	String getVersion();

	boolean hasBeenRun();

	String getSearchType();
}
