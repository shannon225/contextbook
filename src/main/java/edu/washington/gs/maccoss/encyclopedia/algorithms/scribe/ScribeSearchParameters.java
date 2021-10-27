package edu.washington.gs.maccoss.encyclopedia.algorithms.scribe;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import edu.washington.gs.maccoss.encyclopedia.Encyclopedia;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorVersion;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ScoringBreadthType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public class ScribeSearchParameters extends SearchParameters {

	public void savePreferences(File libraryFile, File fastaFile) throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("Scribe");
		HashMap<String, String> map=toParameterMap();
		if (libraryFile!=null) map.put(Encyclopedia.TARGET_LIBRARY_TAG, libraryFile.getAbsolutePath());
		if (fastaFile!=null) map.put(Encyclopedia.BACKGROUND_FASTA_TAG, fastaFile.getAbsolutePath());
		for (Entry<String, String> entry : map.entrySet()) {
			//System.out.println("Writing Scribe preference "+entry.getKey()+" = "+entry.getValue());
			prefs.put(entry.getKey(), entry.getValue());
		}
		prefs.flush();
	}
	
	public static HashMap<String, String> readPreferences() throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("Scribe");
		HashMap<String, String> map=new HashMap<String, String>();
		for (String key : prefs.keys()) {
			String value=prefs.get(key, "");
			//System.out.println("Reading Scribe preference "+key+" = "+value);
			map.put(key, value);
		}
		return map;
	}

	
	public HashMap<String, String> toParameterMap() {
		HashMap<String, String> map=super.toParameterMap();

		return map;
	}
	
	public static ScribeSearchParameters parseParameters(HashMap<String, String> map) {
		SearchParameters params=SearchParameterParser.parseParameters(map);
		return convertFromEncyclopeDIA(params);
	}


	public ScribeSearchParameters(
			AminoAcidConstants aaConstants,
			FragmentationType fragType,
			MassTolerance precursorTolerance,
			double precursorOffsetPPM,
			double precursorIsolationMargin,
			MassTolerance fragmentTolerance,
			double fragmentOffsetPPM,
			MassTolerance libraryFragmentTolerance,
			DigestionEnzyme enzyme,
			float percolatorThreshold,
			float percolatorProteinThreshold,
			PercolatorVersion percolatorVersionNumber,
			int percolatorTrainingSetSize,
			float percolatorTrainingSetThreshold,
			DataAcquisitionType dataAcquisitionType,
			int numberOfThreadsUsed,
			float expectedPeakWidth,
			float targetWindowCenter,
			float precursorWindowSize,
			int numberOfQuantitativePeaks,
			int minNumOfQuantitativePeaks,
			int topNTargetsUsed,
			float minIntensity,
			Optional<PeptideModification> modification,
			ScoringBreadthType searchType,
			float getNumberOfExtraDecoyLibrariesSearched,
			boolean quantifyAcrossSamples,
			boolean verifyModificationIons,
			boolean filterPeaklists,
			boolean doNotUseGlobalFDR,
			Optional<File> precursorIsolationRangeFile, 
			Optional<File> percolatorModelFile,
			boolean enableAdvancedOptions
	) {
		super(
				aaConstants,
				fragType,
				precursorTolerance,
				precursorOffsetPPM,
				precursorIsolationMargin,
				fragmentTolerance,
				fragmentOffsetPPM,
				libraryFragmentTolerance,
				enzyme,
				percolatorThreshold,
				percolatorProteinThreshold,
				percolatorVersionNumber,
				percolatorTrainingSetSize,
				percolatorTrainingSetThreshold,
				dataAcquisitionType,
				numberOfThreadsUsed,
				expectedPeakWidth,
				targetWindowCenter,
				precursorWindowSize,
				numberOfQuantitativePeaks,
				minNumOfQuantitativePeaks,
				topNTargetsUsed,
				minIntensity,
				modification,
				searchType,
				getNumberOfExtraDecoyLibrariesSearched,
				quantifyAcrossSamples,
				verifyModificationIons,
				-1.0f,
				filterPeaklists,
				doNotUseGlobalFDR,
				precursorIsolationRangeFile,
				percolatorModelFile,
				enableAdvancedOptions
		);
	}

	public static ScribeSearchParameters convertFromEncyclopeDIA(SearchParameters params) {
		return new ScribeSearchParameters(
				params.getAAConstants(),
				params.getFragType(),
				params.getPrecursorTolerance(),
				params.getPrecursorOffsetPPM(),
				params.getPrecursorIsolationMargin(),
				params.getFragmentTolerance(),
				params.getFragmentOffsetPPM(),
				params.getLibraryFragmentTolerance(),
				params.getEnzyme(),
				params.getPercolatorThreshold(),
				params.getPercolatorProteinThreshold(),
				params.getPercolatorVersionNumber(),
				params.getPercolatorTrainingSetSize(),
				params.getPercolatorTrainingSetThreshold(),
				params.getDataAcquisitionType(),
				params.getNumberOfThreadsUsed(),
				params.getExpectedPeakWidth(),
				params.getTargetWindowCenter(),
				params.getPrecursorWindowSize(),
				params.getNumberOfQuantitativePeaks(),
				params.getMinNumOfQuantitativePeaks(),
				params.getTopNTargetsUsed(),
				params.getMinIntensity(),
				params.getLocalizingModification(),
				params.getScoringBreadthType(),
				params.getNumberOfExtraDecoyLibrariesSearched(),
				params.isQuantifySameFragmentsAcrossSamples(),
				params.isVerifyModificationIons(),
				params.isFilterPeaklists(),
				params.isDoNotUseGlobalFDR(),
				params.getPrecursorIsolationRangeFile(),
				params.getPercolatorModelFile(),
				params.isEnableAdvancedOptions()
		);
	}
}
