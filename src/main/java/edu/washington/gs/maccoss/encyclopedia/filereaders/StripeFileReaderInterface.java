package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;

public interface StripeFileReaderInterface {
	boolean canTryToReadFile(File file);
    StripeFileInterface createStripeFile(File file, SearchParameters parameters);
}
