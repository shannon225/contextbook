package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import edu.washington.gs.maccoss.encyclopedia.Pecanpie;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public class XCordiaSearchParameters extends PecanSearchParameters {

	public void savePreferences(File backgroundFastaFile, File targetFastaFile) throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("xcordia");
		HashMap<String, String> map=toParameterMap();
		if (backgroundFastaFile!=null) map.put(Pecanpie.BACKGROUND_FASTA_TAG, backgroundFastaFile.getAbsolutePath());
		if (targetFastaFile!=null) map.put(Pecanpie.TARGET_FASTA_TAG, targetFastaFile.getAbsolutePath());
		for (Entry<String, String> entry : map.entrySet()) {
			System.out.println("Writing XCorDIA preference "+entry.getKey()+" = "+entry.getValue());
			prefs.put(entry.getKey(), entry.getValue());
		}
		prefs.flush();
	}
	
	public static HashMap<String, String> readPreferences() throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("xcordia");
		HashMap<String, String> map=new HashMap<String, String>();
		for (String key : prefs.keys()) {
			String value=prefs.get(key, "");
			System.out.println("Reading XCorDIA preference "+key+" = "+value);
			map.put(key, value);
		}
		return map;
	}

	public XCordiaSearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance precursorTolerance, MassTolerance fragmentTolerance, DigestionEnzyme enzyme,
			int maxMissedCleavages, byte minCharge, byte maxCharge, DataAcquisitionType dataAcquisitionType, float precursorWindowSize, int numberOfJobs, int numberOfQuantitativePeaks,
			float numberOfExtraDecoyLibrariesSearched) {
		super(aaConstants, fragType, precursorTolerance, fragmentTolerance, enzyme, maxMissedCleavages, minCharge, maxCharge, dataAcquisitionType, precursorWindowSize, numberOfJobs, numberOfQuantitativePeaks,
				numberOfExtraDecoyLibrariesSearched);
	}

	XCordiaSearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance precursorTolerance, double precursorOffsetPPM, MassTolerance fragmentTolerance,
			double fragmentOffsetPPM, DigestionEnzyme enzyme, int minPeptideLength, int maxPeptideLength, int maxMissedCleavages, byte minCharge, byte maxCharge, int minEluteTime,
			int numberOfReportedPeaks, boolean addDecoysToBackgound, boolean dontRunDecoys, float percolatorThreshold, float alpha, float beta, File percolatorLocation,
			DataAcquisitionType dataAcquisitionType, int numberOfThreadsUsed, float targetWindowCenter, float precursorWindowSize, int numberOfQuantitativePeaks) {
		super(aaConstants, fragType, precursorTolerance, precursorOffsetPPM, fragmentTolerance, fragmentOffsetPPM, enzyme, minPeptideLength, maxPeptideLength, maxMissedCleavages, minCharge, maxCharge,
				minEluteTime, numberOfReportedPeaks, addDecoysToBackgound, dontRunDecoys, percolatorThreshold, alpha, beta, percolatorLocation, dataAcquisitionType, numberOfThreadsUsed, targetWindowCenter,
				precursorWindowSize, numberOfQuantitativePeaks);
	}
	
	public static XCordiaSearchParameters convertFromPecan(PecanSearchParameters params) {
		return new XCordiaSearchParameters(params.getAAConstants(), params.getFragType(), params.getPrecursorTolerance(), params.getFragmentTolerance(), params.getEnzyme(),
				params.getMaxMissedCleavages(), params.getMinCharge(), params.getMaxCharge(), params.getDataAcquisitionType(), params.getPrecursorWindowSize(), 
				params.getNumberOfThreadsUsed(), params.getNumberOfQuantitativePeaks(), params.getNumberOfExtraDecoyLibrariesSearched());
	}

}
