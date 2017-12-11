package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;

public class CASiLJobData extends EncyclopediaJobData {
	public static final String LOG_FILE_SUFFIX=".log";
	public static final String OUTPUT_FILE_SUFFIX=".thesaurus.txt";
	public static final String FEATURE_FILE_SUFFIX=".thesaurus_features.txt";
	public static final String THESAURUS_REPORT_FILE_SUFFIX=".thesaurus" + LibraryFile.ELIB;

	public CASiLJobData(File diaFile, LibraryInterface library, File outputFile, File fastaFile, LibraryScoringFactory taskFactory) {
		super(diaFile, getPercolatorExecutionData(outputFile, fastaFile, taskFactory.getParameters()), taskFactory.getParameters(), taskFactory.getVersion(), library, taskFactory);
	}
	
	private CASiLJobData(File diaFile, PercolatorExecutionData percolatorFiles, SearchParameters parameters, String version, LibraryInterface library, LibraryScoringFactory taskFactory) {
		super(diaFile, percolatorFiles, parameters, version, library, taskFactory);
	}

	static PercolatorExecutionData getPercolatorExecutionData(File referenceFileLocation, File fasta, SearchParameters parameters) {
		return new PercolatorExecutionData(new File(getPrefixFromOutput(referenceFileLocation) + FEATURE_FILE_SUFFIX), fasta,
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_FILE_SUFFIX), 
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_PROTEIN_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_PROTEIN_FILE_SUFFIX), parameters);
	}

	public CASiLJobData updateTaskFactory(LibraryScoringFactory taskFactory) {
		return new CASiLJobData(getDiaFile(), getPercolatorFiles(), getParameters(), getVersion(), getLibrary(), taskFactory);
	}
	
	public File getLocalizationFile() {
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile());
		return new File(absolutePath+".localizations.txt");
	}

	public File getResultLibrary() {
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile());
		return new File(absolutePath+THESAURUS_REPORT_FILE_SUFFIX);
	}

	static String getPrefixFromOutput(File outputFile) {
		final String absolutePath = outputFile.getAbsolutePath();

		if (absolutePath.endsWith(OUTPUT_FILE_SUFFIX)) {
			return absolutePath.substring(0, absolutePath.length() - OUTPUT_FILE_SUFFIX.length());
		} else {
			return absolutePath;
		}
	}
	
	@Override
	public String getSearchType() {
		return "Thesaurus";
	}
}