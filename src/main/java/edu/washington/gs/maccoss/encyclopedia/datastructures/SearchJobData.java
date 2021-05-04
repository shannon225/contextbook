package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;

/**
 * @author Seth Just
 * @since 2:50 PM 9/11/17
 */
public interface SearchJobData {
	String toString();

	StripeFileInterface getDiaFileReader();

	PercolatorExecutionData getPercolatorFiles();

	SearchParameters getParameters();

	String getVersion();

	boolean hasBeenRun();

	String getSearchType();
}
