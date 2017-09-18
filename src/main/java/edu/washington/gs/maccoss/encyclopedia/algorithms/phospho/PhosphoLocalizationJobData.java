package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AbstractSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;

import java.io.File;

public class PhosphoLocalizationJobData extends AbstractSearchJobData {
	public static final String LOG_FILE_SUFFIX=".log";
	public static final String OUTPUT_FILE_SUFFIX=".localization.txt";
	public static final String DECOY_FILE_SUFFIX=".localization.decoy.txt";
	public static final String FEATURE_FILE_SUFFIX=".loc_features.txt";
	
	public PhosphoLocalizationJobData(File diaFile, SearchParameters parameters) {
		super(diaFile, new File(diaFile.getAbsolutePath()+FEATURE_FILE_SUFFIX), new File(diaFile.getAbsolutePath()+OUTPUT_FILE_SUFFIX), new File(diaFile.getAbsolutePath()+DECOY_FILE_SUFFIX), parameters, getJobVersion());
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
		return "PhosphoLocalization";
	}
	public static String getJobVersion() {
		return "0.0.1";
	}
}