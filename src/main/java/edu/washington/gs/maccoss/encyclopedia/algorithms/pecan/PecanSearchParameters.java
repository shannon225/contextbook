package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import edu.washington.gs.maccoss.encyclopedia.Pecanpie;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ScoringBreadthType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

//@Immutable
public class PecanSearchParameters extends SearchParameters {
	private final int minPeptideLength;
	private final int maxPeptideLength;
	private final int maxMissedCleavages;
	private final byte minCharge;
	private final byte maxCharge;
	private final int minEluteTime;
	private final int numberOfReportedPeaks;
	private final boolean addDecoysToBackgound;
	private final boolean dontRunDecoys; // only for testing
	private final float alpha;
	private final float beta;
	public String toString() {
		final StringBuilder sb=new StringBuilder();
		sb.append(" -fixed "+aaConstants.getFixedModString()+"\n");
		sb.append(" -frag "+FragmentationType.toString(fragType)+"\n");
		sb.append(" -ptol "+precursorTolerance.getToleranceThreshold()+"\n");
		sb.append(" -ftol "+fragmentTolerance.getToleranceThreshold()+"\n");
		sb.append(" -ptolunits"+precursorTolerance.getUnits()+"\n");
		sb.append(" -ftolunits"+fragmentTolerance.getUnits()+"\n");
		sb.append(" -poffset "+precursorOffsetPPM+"\n");
		sb.append(" -foffset "+fragmentOffsetPPM+"\n");
		sb.append(" -enzyme "+enzyme.getName()+"\n");
		sb.append(" -minLength "+minPeptideLength+"\n");
		sb.append(" -maxLength "+maxPeptideLength+"\n");
		sb.append(" -maxMissedCleavage "+maxMissedCleavages+"\n");
		sb.append(" -minCharge "+minCharge+"\n");
		sb.append(" -maxCharge "+maxCharge+"\n");
		sb.append(" -minEluteTime "+minEluteTime+"\n");
		sb.append(" -numberOfReportedPeaks "+numberOfReportedPeaks+"\n");
		sb.append(" -addDecoysToBackground "+addDecoysToBackgound+"\n");
		sb.append(" -dontRunDecoys "+dontRunDecoys+"\n");
		sb.append(" -percolatorThreshold "+percolatorThreshold+"\n");
		sb.append(" -percolatorVersionNumber "+percolatorVersionNumber+"\n");
		sb.append(" -acquisition "+DataAcquisitionType.toString(dataAcquisitionType)+"\n");
		sb.append(" -numberOfThreadsUsed "+numberOfThreadsUsed+"\n");
		sb.append(" -numberOfQuantitativePeaks "+numberOfQuantitativePeaks+"\n");
		sb.append(" -minNumOfQuantitativePeaks "+minNumOfQuantitativePeaks+"\n");
		sb.append(" -alpha "+alpha+"\n");
		sb.append(" -beta "+beta+"\n");
		
		return sb.toString();
	}
	
	public HashMap<String, String> toParameterMap() {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("-fixed", aaConstants.getFixedModString());
		map.put("-frag", FragmentationType.toString(fragType));
		map.put("-ptol", precursorTolerance.getToleranceThreshold()+"");
		map.put("-ftol", fragmentTolerance.getToleranceThreshold()+"");
		map.put("-ptolunits", precursorTolerance.getUnits());
		map.put("-ftolunits", fragmentTolerance.getUnits());
		map.put("-poffset", precursorOffsetPPM+"");
		map.put("-foffset", fragmentOffsetPPM+"");
		map.put("-enzyme", enzyme.getName());
		map.put("-minLength", minPeptideLength+"");
		map.put("-maxLength", maxPeptideLength+"");
		map.put("-maxMissedCleavage", maxMissedCleavages+"");
		map.put("-minCharge", minCharge+"");
		map.put("-maxCharge", maxCharge+"");
		map.put("-minEluteTime", minEluteTime+"");
		map.put("-numberOfReportedPeaks", numberOfReportedPeaks+"");
		map.put("-addDecoysToBackground", addDecoysToBackgound+"");
		map.put("-dontRunDecoys", dontRunDecoys+"");
		map.put("-percolatorThreshold", percolatorThreshold+"");
		map.put("-percolatorVersionNumber", percolatorVersionNumber+"");
		map.put("-acquisition", DataAcquisitionType.toString(dataAcquisitionType));
		map.put("-numberOfThreadsUsed", numberOfThreadsUsed+"");
		map.put("-precursorWindowSize", precursorWindowSize+"");
		map.put("-numberOfQuantitativePeaks", numberOfQuantitativePeaks+"");
		map.put("-minNumOfQuantitativePeaks", minNumOfQuantitativePeaks+"");
		map.put("-alpha", alpha+"");
		map.put("-beta", beta+"");
		return map;
	}
	
	public void savePreferences(File backgroundFastaFile, File targetFastaFile) throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("pecan");
		HashMap<String, String> map=toParameterMap();
		if (backgroundFastaFile!=null) map.put(Pecanpie.BACKGROUND_FASTA_TAG, backgroundFastaFile.getAbsolutePath());
		if (targetFastaFile!=null) map.put(Pecanpie.TARGET_FASTA_TAG, targetFastaFile.getAbsolutePath());
		for (Entry<String, String> entry : map.entrySet()) {
			//System.out.println("Writing Pecan preference "+entry.getKey()+" = "+entry.getValue());
			prefs.put(entry.getKey(), entry.getValue());
		}
		prefs.flush();
	}
	
	public static HashMap<String, String> readPreferences() throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("pecan");
		HashMap<String, String> map=new HashMap<String, String>();
		for (String key : prefs.keys()) {
			String value=prefs.get(key, "");
			//System.out.println("Reading Pecan preference "+key+" = "+value);
			map.put(key, value);
		}
		return map;
	}
	
	public PecanSearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance precursorTolerance, double precursorOffsetPPM, double precursorIsolationMargin, MassTolerance fragmentTolerance,
			double fragmentOffsetPPM, DigestionEnzyme enzyme, int minPeptideLength, int maxPeptideLength, int maxMissedCleavages, byte minCharge, byte maxCharge, int minEluteTime,
			int numberOfReportedPeaks, boolean addDecoysToBackgound, boolean dontRunDecoys, float percolatorThreshold, float percolatorProteinThreshold, float alpha, float beta, Integer percolatorVersionNumber,
			DataAcquisitionType dataAcquisitionType, int numberOfThreadsUsed, float targetWindowCenter, float precursorWindowSize, int numberOfQuantitativePeaks, int minNumOfQuantitativePeaks, int minQuantitativeIonNumber, boolean quantifyAcrossSamples) {
		super(aaConstants, fragType, precursorTolerance, precursorOffsetPPM, precursorIsolationMargin, fragmentTolerance, fragmentOffsetPPM, fragmentTolerance, enzyme, percolatorThreshold, percolatorProteinThreshold, percolatorVersionNumber, dataAcquisitionType, numberOfThreadsUsed, minEluteTime*2.0f, targetWindowCenter, precursorWindowSize, numberOfQuantitativePeaks, minNumOfQuantitativePeaks, minQuantitativeIonNumber, Optional.ofNullable((PeptideModification)null), ScoringBreadthType.ENTIRE_RT_WINDOW, 0, quantifyAcrossSamples);
		this.minPeptideLength=minPeptideLength;
		this.maxPeptideLength=maxPeptideLength;
		this.maxMissedCleavages=maxMissedCleavages;
		this.minCharge=minCharge;
		this.maxCharge=maxCharge;
		this.minEluteTime=minEluteTime;
		this.numberOfReportedPeaks=numberOfReportedPeaks;
		this.addDecoysToBackgound=addDecoysToBackgound;
		this.dontRunDecoys=dontRunDecoys;
		this.alpha=alpha;
		this.beta=beta;
	}
	
	public PecanSearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance precursorTolerance, MassTolerance fragmentTolerance, DigestionEnzyme enzyme, int percolatorVersionNumber,
			int maxMissedCleavages, byte minCharge, byte maxCharge, DataAcquisitionType dataAcquisitionType, float precursorWindowSize, int numberOfJobs, int numberOfQuantitativePeaks, int minNumOfQuantitativePeaks, int minQuantitativeIonNumber, float numberOfExtraDecoyLibrariesSearched, boolean quantifyAcrossSamples) {
		super(aaConstants, fragType, precursorTolerance, 0.0, 0.0, fragmentTolerance, 0.0, fragmentTolerance, enzyme, 0.01f, 0.01f, percolatorVersionNumber, dataAcquisitionType, numberOfJobs, 24f, -1f, precursorWindowSize, numberOfQuantitativePeaks, minNumOfQuantitativePeaks, minQuantitativeIonNumber, Optional.ofNullable((PeptideModification)null), ScoringBreadthType.ENTIRE_RT_WINDOW, numberOfExtraDecoyLibrariesSearched, quantifyAcrossSamples);
		minPeptideLength=5;
		maxPeptideLength=100;
		this.maxMissedCleavages=maxMissedCleavages;
		this.minCharge=minCharge;
		this.maxCharge=maxCharge;
		minEluteTime=12;
		numberOfReportedPeaks=1;
		addDecoysToBackgound=false;
		dontRunDecoys=false;
		alpha=1.8f;
		beta=0.4f;
	}

	public PecanSearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance fragmentTolerance, MassTolerance precursorTolerance, DigestionEnzyme enzyme, DataAcquisitionType dataAcquisitionType, boolean quantifyAcrossSamples) {
		super(aaConstants, fragType, precursorTolerance, 0.0, 0.0, fragmentTolerance, 0.0, fragmentTolerance, enzyme, 0.01f, 0.01f, null, dataAcquisitionType, Runtime.getRuntime().availableProcessors(), 24f, -1f, -1f, 5, 0, 0, Optional.ofNullable((PeptideModification)null), ScoringBreadthType.ENTIRE_RT_WINDOW, 0, quantifyAcrossSamples);
		minPeptideLength=5;
		maxPeptideLength=100;
		maxMissedCleavages=1;
		minCharge=2;
		maxCharge=3;
		minEluteTime=12;
		numberOfReportedPeaks=1;
		addDecoysToBackgound=false;
		dontRunDecoys=false;
		alpha=1.8f;
		beta=0.4f;
	}

	public PecanSearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance fragmentTolerance, MassTolerance precursorTolerance, DigestionEnzyme enzyme, boolean quantifyAcrossSamples) {
		super(aaConstants, fragType, precursorTolerance, 0.0, 0.0, fragmentTolerance, 0.0, fragmentTolerance, enzyme, 0.01f, 0.01f, null, DataAcquisitionType.DIA, Runtime.getRuntime().availableProcessors(), 24f, -1f, -1f, 5, 0, 0, Optional.ofNullable((PeptideModification)null), ScoringBreadthType.ENTIRE_RT_WINDOW, 0, quantifyAcrossSamples);
		minPeptideLength=5;
		maxPeptideLength=100;
		maxMissedCleavages=1;
		minCharge=2;
		maxCharge=3;
		minEluteTime=12;
		numberOfReportedPeaks=1;
		addDecoysToBackgound=false;
		dontRunDecoys=false;
		alpha=1.8f;
		beta=0.4f;
	}

	public PecanSearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance fragmentTolerance, MassTolerance precursorTolerance, DigestionEnzyme enzyme,
			int maxMissedCleavages, boolean quantifyAcrossSamples) {
		super(aaConstants, fragType, precursorTolerance, 0.0, 0.0, fragmentTolerance, 0.0, fragmentTolerance, enzyme, 0.01f, 0.01f, null, DataAcquisitionType.DIA, Runtime.getRuntime().availableProcessors(), 24f, -1f, -1f, 5, 0, 0, Optional.ofNullable((PeptideModification)null), ScoringBreadthType.ENTIRE_RT_WINDOW, 0, quantifyAcrossSamples);
		this.maxMissedCleavages=maxMissedCleavages;
		minPeptideLength=5;
		maxPeptideLength=100;
		minCharge=2;
		maxCharge=3;
		minEluteTime=12;
		numberOfReportedPeaks=1;
		addDecoysToBackgound=false;
		dontRunDecoys=false;
		alpha=1.8f;
		beta=0.4f;
	}

	public int getMaxMissedCleavages() {
		return maxMissedCleavages;
	}

	public int getMaxPeptideLength() {
		return maxPeptideLength;
	}

	public int getMinPeptideLength() {
		return minPeptideLength;
	}

	public byte getMaxCharge() {
		return maxCharge;
	}

	public byte getMinCharge() {
		return minCharge;
	}

	public int getMinEluteTime() {
		return minEluteTime;
	}

	public int getNumberOfReportedPeaks() {
		return numberOfReportedPeaks;
	}

	public boolean isAddDecoysToBackgound() {
		return addDecoysToBackgound;
	}
	
	public boolean isDontRunDecoys() {
		return dontRunDecoys;
	}
	
	public float getAlpha() {
		return alpha;
	}
	
	public float getBeta() {
		return beta;
	}
}
