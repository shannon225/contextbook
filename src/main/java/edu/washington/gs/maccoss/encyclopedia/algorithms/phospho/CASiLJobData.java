package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;

public class CASiLJobData extends EncyclopediaJobData {
	public static final String LOG_FILE_SUFFIX=".log";
	public static final String OUTPUT_FILE_SUFFIX=".casil.txt";
	public static final String FEATURE_FILE_SUFFIX=".casil_features.txt";
	
	public CASiLJobData(File diaFile, LibraryInterface library, LibraryScoringFactory taskFactory) {
		super(diaFile, new File(diaFile.getAbsolutePath()+FEATURE_FILE_SUFFIX), new File(diaFile.getAbsolutePath()+OUTPUT_FILE_SUFFIX), taskFactory.getParameters(), taskFactory.getVersion(), library, taskFactory);
	}

	public CASiLJobData(File diaFile, LibraryInterface library, File outputFile, LibraryScoringFactory taskFactory) {
		super(diaFile, new File(getOutputAbsolutePathPrefix(outputFile.getAbsolutePath())+FEATURE_FILE_SUFFIX), outputFile, taskFactory.getParameters(), taskFactory.getVersion(), library, taskFactory);
	}
	
	private CASiLJobData(File diaFile, File featureFile, File outputFile, SearchParameters parameters, String version, LibraryInterface library, LibraryScoringFactory taskFactory) {
		super(diaFile, featureFile, outputFile, parameters, version, library, taskFactory);
	}

	public CASiLJobData updateTaskFactory(LibraryScoringFactory taskFactory) {
		return new CASiLJobData(getDiaFile(), getFeatureFile(), getOutputFile(), getParameters(), getVersion(), getLibrary(), taskFactory);
	}
	
	public File getFirstPassPercolator() {
		// note, overwrites the first pass with the actual output
		return getOutputFile();
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
		return "CASiL";
	}
}