package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class ThesaurusJobData extends EncyclopediaJobData {
	public static final String LOG_FILE_SUFFIX=".log";
	public static final String OUTPUT_FILE_SUFFIX=".thesaurus.txt";
	public static final String FEATURE_FILE_SUFFIX=".thesaurus_features.txt";
	public static final String THESAURUS_REPORT_FILE_SUFFIX=".thesaurus" + LibraryFile.ELIB;

	public ThesaurusJobData(File diaFile, LibraryInterface library, File outputFile, File fastaFile, LibraryScoringFactory taskFactory) {
		super(diaFile, getPercolatorExecutionData(outputFile, fastaFile, taskFactory.getParameters()), taskFactory.getParameters(), taskFactory.getVersion(), library, taskFactory);
	}
	
	private ThesaurusJobData(File diaFile, PercolatorExecutionData percolatorFiles, SearchParameters parameters, String version, LibraryInterface library, LibraryScoringFactory taskFactory) {
		super(diaFile, percolatorFiles, parameters, version, library, taskFactory);
	}

	static PercolatorExecutionData getPercolatorExecutionData(File referenceFileLocation, File fasta, SearchParameters parameters) {
		String prefix = getPrefix(parameters);
		return new PercolatorExecutionData(new File(getPrefixFromOutput(referenceFileLocation, parameters) + prefix+FEATURE_FILE_SUFFIX), fasta,
				new File(getPrefixFromOutput(referenceFileLocation, parameters) + prefix+OUTPUT_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation, parameters) + DECOY_FILE_SUFFIX), 
				new File(getPrefixFromOutput(referenceFileLocation, parameters) + OUTPUT_PROTEIN_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation, parameters) + DECOY_PROTEIN_FILE_SUFFIX), parameters);
	}

	public ThesaurusJobData updateTaskFactory(LibraryScoringFactory taskFactory) {
		return new ThesaurusJobData(getDiaFile(), getPercolatorFiles(), getParameters(), getVersion(), getLibrary(), taskFactory);
	}
	
	public File getLocalizationFile() {
		String prefix = getPrefix(getParameters());
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile(), getParameters());
		return new File(absolutePath+prefix+".localizations.txt");
	}

	public File getResultLibrary() {
		String prefix = getPrefix(getParameters());
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile(), getParameters());
		return new File(absolutePath+prefix+THESAURUS_REPORT_FILE_SUFFIX);
	}

	static String getPrefixFromOutput(File outputFile, SearchParameters parameters) {
		final String absolutePath = outputFile.getAbsolutePath();

		String prefix = getPrefix(parameters);
		if (absolutePath.endsWith(prefix+OUTPUT_FILE_SUFFIX)) {
			return absolutePath.substring(0, absolutePath.length() - (prefix.length()+OUTPUT_FILE_SUFFIX.length()));
		} else {
			return absolutePath;
		}
	}
	
	@Override
	public String getSearchType() {
		return "Thesaurus";
	}

	private static String getPrefix(SearchParameters parameters) {
		if (!parameters.getLocalizingModification().isPresent()) {
			throw new EncyclopediaException("Sorry, Thesaurus requires you to specify a localizing PTM!");
		}
		String prefix="."+parameters.getLocalizingModification().get().getShortname();
		return prefix;
	}
}