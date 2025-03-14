package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Optional;
import java.util.Map.Entry;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import edu.washington.gs.maccoss.encyclopedia.Pecanpie;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorVersion;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.parameters.InstrumentSpecificSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public class XCordiaSearchParameters extends PecanSearchParameters {
	@Override
	public void savePreferences(File backgroundFastaFile, File targetFastaFile) throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("xcordia");
		HashMap<String, String> map=toParameterMap();
		if (backgroundFastaFile!=null) map.put(Pecanpie.BACKGROUND_FASTA_TAG, backgroundFastaFile.getAbsolutePath());
		if (targetFastaFile!=null) map.put(Pecanpie.TARGET_FASTA_TAG, targetFastaFile.getAbsolutePath());
		for (Entry<String, String> entry : map.entrySet()) {
			//System.out.println("Writing XCorDIA preference "+entry.getKey()+" = "+entry.getValue());
			prefs.put(entry.getKey(), entry.getValue());
		}
		prefs.flush();
	}
	
	public static HashMap<String, String> readPreferences() throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("xcordia");
		HashMap<String, String> map=new HashMap<String, String>();
		for (String key : prefs.keys()) {
			String value=prefs.get(key, "");
			//System.out.println("Reading XCorDIA preference "+key+" = "+value);
			map.put(key, value);
		}
		return map;
	}

	/** used by CLI
	 */
	public XCordiaSearchParameters(
			AminoAcidConstants aaConstants,
			FragmentationType fragType,
			MassTolerance precursorTolerance,
			double precursorOffsetPPM,
			double precursorIsolationMargin,
			MassTolerance fragmentTolerance,
			double fragmentOffsetPPM,
			DigestionEnzyme enzyme,
			float expectedPeakWidth,
			int minPeptideLength,
			int maxPeptideLength,
			int maxMissedCleavages,
			byte minCharge,
			byte maxCharge,
			int numberOfReportedPeaks,
			boolean addDecoysToBackgound,
			boolean dontRunDecoys,
			float percolatorThreshold,
			float percolatorProteinThreshold,
			boolean usePercolator,
			PercolatorVersion percolatorVersionNumber,
			int percolatorTrainingSetSize,
			float percolatorTrainingSetThreshold,
			int percolatorTrainingIterations,
			float alpha,
			float beta,
			DataAcquisitionType dataAcquisitionType,
			int numberOfThreadsUsed,
			float targetWindowCenter,
			float precursorWindowSize,
			float maxWindowWidth,
			int numberOfQuantitativePeaks,
			int minNumOfQuantitativePeaks,
			int topNTargetsUsed,
			float minIntensity,
			float minIntensityNumIons, 
			boolean quantifyAcrossSamples,
			boolean verifyModificationIons,
			boolean requireVariableMods,
			float rtWindowInMin,
			int minNumIntegratedRTPoints,
			boolean filterPeaklists,
			boolean doNotUseGlobalFDR,
			Optional<File> precursorIsolationRangeFile, 
			Optional<File> percolatorModelFile, 
			boolean smoothIntegrations,
			boolean normalizeByTIC,
			boolean subtractBackground,
			boolean maskBadIntegrations,
			boolean adjustInferredRTBoundaries,
			boolean skipLibraryRetentionTime,
			boolean integratePrecursors,
			InstrumentSpecificSearchParameters instrument,
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
				enzyme,
				expectedPeakWidth,
				minPeptideLength,
				maxPeptideLength,
				maxMissedCleavages,
				minCharge,
				maxCharge,
				numberOfReportedPeaks,
				addDecoysToBackgound,
				dontRunDecoys,
				percolatorThreshold,
				percolatorProteinThreshold,
				usePercolator,
				percolatorVersionNumber,
				percolatorTrainingSetSize,
				percolatorTrainingSetThreshold,
				percolatorTrainingIterations,
				alpha,
				beta,
				dataAcquisitionType,
				numberOfThreadsUsed,
				targetWindowCenter,
				precursorWindowSize,
				maxWindowWidth,
				numberOfQuantitativePeaks,
				minNumOfQuantitativePeaks,
				topNTargetsUsed,
				minIntensity,
				minIntensityNumIons, 
				quantifyAcrossSamples,
				verifyModificationIons,
				requireVariableMods,
				rtWindowInMin,
				minNumIntegratedRTPoints,
				filterPeaklists,
				doNotUseGlobalFDR,
				precursorIsolationRangeFile,
				percolatorModelFile,
				smoothIntegrations,
				normalizeByTIC,
				subtractBackground,
				maskBadIntegrations,
				adjustInferredRTBoundaries,
				skipLibraryRetentionTime,
				integratePrecursors,
				instrument,
				enableAdvancedOptions
		);
	}

	/** used by GUI
	 */
	public XCordiaSearchParameters(
			AminoAcidConstants aaConstants,
			FragmentationType fragType,
			MassTolerance precursorTolerance,
			MassTolerance fragmentTolerance,
			DigestionEnzyme enzyme,
			PercolatorVersion percolatorVersionNumber,
			float percolatorThreshold,
			float percolatorProteinThreshold,
			int percolatorTrainingSetSize,
			float percolatorTrainingSetThreshold,
			int percolatorTrainingIterations,
			int maxMissedCleavages,
			byte minCharge,
			byte maxCharge,
			DataAcquisitionType dataAcquisitionType,
			float precursorWindowSize,
			int numberOfJobs,
			int numberOfQuantitativePeaks,
			int minNumOfQuantitativeIons,
			int topNTargetsUsed,
			float minIntensity,
			float numberOfExtraDecoyLibrariesSearched,
			boolean quantifyAcrossSamples,
			boolean verifyModificationIons,
			boolean requireVariableMods,
			InstrumentSpecificSearchParameters instrument
	) {
		super(
				aaConstants,
				fragType,
				precursorTolerance,
				fragmentTolerance,
				enzyme,
				percolatorVersionNumber,
				percolatorThreshold,
				percolatorProteinThreshold,
				percolatorTrainingSetSize,
				percolatorTrainingSetThreshold,
				percolatorTrainingIterations,
				maxMissedCleavages,
				minCharge,
				maxCharge,
				dataAcquisitionType,
				precursorWindowSize,
				numberOfJobs,
				numberOfQuantitativePeaks,
				minNumOfQuantitativeIons,
				topNTargetsUsed,
				minIntensity,
				numberOfExtraDecoyLibrariesSearched,
				quantifyAcrossSamples,
				verifyModificationIons,
				requireVariableMods,
				instrument
		);
	}
	
	public static XCordiaSearchParameters convertFromPecan(PecanSearchParameters params) {
		return new XCordiaSearchParameters(
				params.getAAConstants(),
				params.getFragType(),
				params.getPrecursorTolerance(),
				params.getPrecursorOffsetPPM(),
				params.getPrecursorIsolationMargin(),
				params.getFragmentTolerance(),
				params.getFragmentOffsetPPM(),
				params.getEnzyme(),
				params.getExpectedPeakWidth(),
				params.getMinPeptideLength(),
				params.getMaxPeptideLength(),
				params.getMaxMissedCleavages(),
				params.getMinCharge(),
				params.getMaxCharge(),
				params.getNumberOfReportedPeaks(),
				params.isAddDecoysToBackgound(),
				params.isDontRunDecoys(),
				params.getPercolatorThreshold(),
				params.getPercolatorProteinThreshold(),
				params.isUsePercolator(),
				params.getPercolatorVersionNumber(),
				params.getPercolatorTrainingSetSize(),
				params.getPercolatorTrainingSetThreshold(),
				params.getPercolatorTrainingIterations(),
				params.getAlpha(),
				params.getBeta(),
				params.getDataAcquisitionType(),
				params.getNumberOfThreadsUsed(),
				params.getTargetWindowCenter(),
				params.getPrecursorWindowSize(),
				params.getMaxWindowWidth(),
				params.getNumberOfQuantitativePeaks(),
				params.getMinNumOfQuantitativePeaks(),
				params.getTopNTargetsUsed(),
				params.getMinIntensity(),
				params.getMinIntensityNumIons(),
				params.isQuantifySameFragmentsAcrossSamples(),
				params.isVerifyModificationIons(),
				params.isRequireVariableMods(),
				params.getRtWindowInMin(),
				params.getMinNumIntegratedRTPoints(),
				params.isFilterPeaklists(),
				params.isDoNotUseGlobalFDR(),
				params.getPrecursorIsolationRangeFile(),
				params.getPercolatorModelFile(),
				params.isSmoothIntegrations(),
				params.isNormalizeByTIC(),
				params.isSubtractBackground(),
				params.isMaskBadIntegrations(),
				params.adjustInferredRTBoundaries(), 
				params.isSkipLibraryRetentionTime(),
				params.isIntegratePrecursors(),
				params.getInstrument(),
				params.isEnableAdvancedOptions()
		);
	}

}
