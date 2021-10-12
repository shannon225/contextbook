package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import edu.washington.gs.maccoss.encyclopedia.Encyclopedia;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorVersion;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public class ThesaurusSearchParameters extends SearchParameters {
	private final boolean considerRearrangement; // default is false

	public void savePreferences(File libraryFile, File fastaFile) throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("CASiL");
		HashMap<String, String> map=toParameterMap();
		if (libraryFile!=null) map.put(Encyclopedia.TARGET_LIBRARY_TAG, libraryFile.getAbsolutePath());
		if (fastaFile!=null) map.put(Encyclopedia.BACKGROUND_FASTA_TAG, fastaFile.getAbsolutePath());
		for (Entry<String, String> entry : map.entrySet()) {
			//System.out.println("Writing EncyclopeDIA preference "+entry.getKey()+" = "+entry.getValue());
			prefs.put(entry.getKey(), entry.getValue());
		}
		prefs.flush();
	}
	
	public static HashMap<String, String> readPreferences() throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("CASiL");
		HashMap<String, String> map=new HashMap<String, String>();
		for (String key : prefs.keys()) {
			String value=prefs.get(key, "");
			//System.out.println("Reading CASiL preference "+key+" = "+value);
			map.put(key, value);
		}
		return map;
	}

	
	public HashMap<String, String> toParameterMap() {
		HashMap<String, String> map=super.toParameterMap();

        map.put("-considerRearrangement", considerRearrangement+"");
		return map;
	}
	
	public static ThesaurusSearchParameters parseParameters(HashMap<String, String> map) {
		SearchParameters params=SearchParameterParser.parseParameters(map);
		String value=map.get("-considerRearrangement");
		boolean considerRearrangement=value!=null&&Boolean.getBoolean(value); // only true if not null and "true"
		return convertFromEncyclopeDIA(params, considerRearrangement);
	}


	public ThesaurusSearchParameters(
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
			PeptideModification modification,
			ScoringBreadthType searchType,
			float getNumberOfExtraDecoyLibrariesSearched,
			boolean quantifyAcrossSamples,
			boolean verifyModificationIons,
			boolean filterPeaklists,
			boolean doNotUseGlobalFDR,
			boolean considerRearrangement, 
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
				Optional.of(modification),
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
		this.considerRearrangement=considerRearrangement;
	}

	public static ThesaurusSearchParameters convertFromEncyclopeDIA(SearchParameters params, boolean considerRearrangement) {
		PeptideModification mod;
		if (params.getLocalizingModification().isPresent()) {
			mod=params.getLocalizingModification().get();
		} else {
			Logger.logLine("You should specify a localization modification if you're going to apply localization! Using phosphorylation by default.");
			mod=PeptideModification.phosphorylation;
		}
		return new ThesaurusSearchParameters(
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
				mod,
				params.getScoringBreadthType(),
				params.getNumberOfExtraDecoyLibrariesSearched(),
				params.isQuantifySameFragmentsAcrossSamples(),
				params.isVerifyModificationIons(),
				params.isFilterPeaklists(),
				params.isDoNotUseGlobalFDR(),
				considerRearrangement,
				params.getPrecursorIsolationRangeFile(),
				params.getPercolatorModelFile(),
				params.isEnableAdvancedOptions()
		);
	}

	public boolean isConsiderRearrangement() {
		return considerRearrangement;
	}
}
