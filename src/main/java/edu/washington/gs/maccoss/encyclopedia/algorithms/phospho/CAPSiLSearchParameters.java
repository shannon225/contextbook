package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import edu.washington.gs.maccoss.encyclopedia.Encyclopedia;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public class CAPSiLSearchParameters extends SearchParameters {

	
	public void savePreferences(File libraryFile) throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("capsil");
		HashMap<String, String> map=toParameterMap();
		if (libraryFile!=null) map.put(Encyclopedia.TARGET_LIBRARY_TAG, libraryFile.getAbsolutePath());
		for (Entry<String, String> entry : map.entrySet()) {
			//System.out.println("Writing EncyclopeDIA preference "+entry.getKey()+" = "+entry.getValue());
			prefs.put(entry.getKey(), entry.getValue());
		}
		prefs.flush();
	}
	
	public static HashMap<String, String> readPreferences() throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("capsil");
		HashMap<String, String> map=new HashMap<String, String>();
		for (String key : prefs.keys()) {
			String value=prefs.get(key, "");
			//System.out.println("Reading EncyclopeDIA preference "+key+" = "+value);
			map.put(key, value);
		}
		return map;
	}
	
	

	public CAPSiLSearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance precursorTolerance, double precursorOffsetPPM, double precursorIsolationMargin,
			MassTolerance fragmentTolerance, double fragmentOffsetPPM, MassTolerance libraryFragmentTolerance, DigestionEnzyme enzyme, float percolatorThreshold, Integer percolatorVersionNumber,
			DataAcquisitionType dataAcquisitionType, int numberOfThreadsUsed, float expectedPeakWidth, float targetWindowCenter, float precursorWindowSize, int numberOfQuantitativePeaks,
			int minNumOfQuantitativePeaks, PeptideModification modification, ScoringBreadthType searchType, float getNumberOfExtraDecoyLibrariesSearched) {
		super(aaConstants, fragType, precursorTolerance, precursorOffsetPPM, precursorIsolationMargin, fragmentTolerance, fragmentOffsetPPM, libraryFragmentTolerance, enzyme, percolatorThreshold,
				percolatorVersionNumber, dataAcquisitionType, numberOfThreadsUsed, expectedPeakWidth, targetWindowCenter, precursorWindowSize, numberOfQuantitativePeaks, minNumOfQuantitativePeaks,
				Optional.of(modification), searchType, getNumberOfExtraDecoyLibrariesSearched);
	}

	public static CAPSiLSearchParameters convertFromEncyclopeDIA(SearchParameters params) {
		PeptideModification mod;
		if (params.getLocalizingModification().isPresent()) {
			mod=params.getLocalizingModification().get();
		} else {
			Logger.logLine("You should specify a localization modification if you're going to apply localization! Using phosphorylation by default.");
			mod=PeptideModification.phosphorylation;
		}
		return new CAPSiLSearchParameters(params.getAAConstants(), params.getFragType(), params.getPrecursorTolerance(), params.getPrecursorOffsetPPM(), params.getPrecursorIsolationMargin(),
				params.getFragmentTolerance(), params.getFragmentOffsetPPM(), params.getLibraryFragmentTolerance(), params.getEnzyme(), params.getPercolatorThreshold(),
				params.getPercolatorVersionNumber(), params.getDataAcquisitionType(), params.getNumberOfThreadsUsed(), params.getExpectedPeakWidth(), params.getTargetWindowCenter(),
				params.getPrecursorWindowSize(), params.getNumberOfQuantitativePeaks(), params.getMinNumOfQuantitativePeaks(), mod, params.getScoringBreadthType(),
				params.getNumberOfExtraDecoyLibrariesSearched());
	}

}
