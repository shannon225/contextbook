package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.FilenameFilter;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;

public interface StripeFileReaderInterface extends FilenameFilter {
	boolean canTryToReadFile(File file);
    StripeFileInterface readStripeFile(File file, SearchParameters parameters, boolean isOpenFileInPlace);
}
