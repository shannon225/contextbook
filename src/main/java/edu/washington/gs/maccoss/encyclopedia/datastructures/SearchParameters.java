package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import edu.washington.gs.maccoss.encyclopedia.Encyclopedia;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ScoringBreadthType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public class SearchParameters {
	public static final String OPT_PERC_TRAINING_SIZE = "-percolatorTrainingSetSize";
	public static final String OPT_PERC_TRAINING_THRESH = "-percolatorTrainingFDR";

	protected final AminoAcidConstants aaConstants;
	protected final FragmentationType fragType;
	protected final MassTolerance precursorTolerance;
	protected final MassTolerance fragmentTolerance;
	protected final MassTolerance libraryFragmentTolerance;
	protected final DigestionEnzyme enzyme;
	protected final float percolatorThreshold;
	protected final float percolatorProteinThreshold;
	protected final int percolatorVersionNumber;
	protected final int percolatorTrainingSetSize;
	protected final float percolatorTrainingSetThreshold;
	protected final DataAcquisitionType dataAcquisitionType;
	protected final int numberOfThreadsUsed;	
	protected final float targetWindowCenter;
	protected final float expectedPeakWidth;
	protected final float precursorWindowSize;
	protected final float numberOfExtraDecoyLibrariesSearched;
	protected final int numberOfQuantitativePeaks;
	protected final int minNumOfQuantitativePeaks;
	protected final float minIntensity;
	protected final double precursorOffsetPPM;
	protected final double fragmentOffsetPPM;
	protected final double precursorIsolationMargin;
	protected final boolean useNLsForXCorr=false;
	protected final ScoringBreadthType CASiLBreadthType;
	protected final Optional<PeptideModification> localizingModification; 
	protected final boolean quantifyAcrossSamples;
	protected final boolean verifyModificationIons;
    protected final float rtWindowInMin;
    protected final boolean filterPeaklists;
    protected final boolean doNotUseGlobalFDR;
    protected final int topNTargetsUsed;

	public SearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance precursorTolerance, double precursorOffsetPPM, double precursorIsolationMargin, MassTolerance fragmentTolerance, double fragmentOffsetPPM, MassTolerance libraryFragmentTolerance, DigestionEnzyme enzyme,
			float percolatorThreshold, float percolatorProteinThreshold, Integer percolatorVersionNumber, int percolatorTrainingSetSize, float percolatorTrainingSetThreshold,
			DataAcquisitionType dataAcquisitionType, int numberOfThreadsUsed, float expectedPeakWidth, float targetWindowCenter, float precursorWindowSize,
			int numberOfQuantitativePeaks, int minNumOfQuantitativePeaks, int topNTargetsUsed, float minIntensity, Optional<PeptideModification> localizingModification, ScoringBreadthType CASiLBreadthType, float getNumberOfExtraDecoyLibrariesSearched, boolean quantifyAcrossSamples, boolean verifyModificationIons, float rtWindowInMin, boolean filterPeaklists, boolean doNotUseGlobalFDR) {
		this.aaConstants=aaConstants;
		this.fragType=fragType;
		this.precursorTolerance=precursorTolerance;
		this.precursorOffsetPPM=precursorOffsetPPM;
		this.precursorIsolationMargin=precursorIsolationMargin;
		this.fragmentTolerance=fragmentTolerance;
		this.fragmentOffsetPPM=fragmentOffsetPPM;
		this.libraryFragmentTolerance=libraryFragmentTolerance;
		this.enzyme=enzyme;
		this.percolatorThreshold=percolatorThreshold;
		this.percolatorProteinThreshold=percolatorProteinThreshold;
		this.percolatorVersionNumber=percolatorVersionNumber==null?PercolatorExecutor.DEFAULT_VERSION_NUMBER:percolatorVersionNumber;
		this.percolatorTrainingSetSize = percolatorTrainingSetSize;
		this.percolatorTrainingSetThreshold = percolatorTrainingSetThreshold;
		this.dataAcquisitionType=dataAcquisitionType;
		this.numberOfThreadsUsed=numberOfThreadsUsed;
		this.expectedPeakWidth=expectedPeakWidth;
		this.targetWindowCenter=targetWindowCenter;
		this.precursorWindowSize=precursorWindowSize;
		this.numberOfQuantitativePeaks=numberOfQuantitativePeaks;
		this.minNumOfQuantitativePeaks=minNumOfQuantitativePeaks;
		this.topNTargetsUsed=topNTargetsUsed;
		this.minIntensity=minIntensity;
		this.CASiLBreadthType=CASiLBreadthType;
		this.localizingModification=localizingModification;
		this.numberOfExtraDecoyLibrariesSearched=getNumberOfExtraDecoyLibrariesSearched;
		this.quantifyAcrossSamples=quantifyAcrossSamples;
		this.verifyModificationIons=verifyModificationIons;
        this.rtWindowInMin=rtWindowInMin;
        this.filterPeaklists=filterPeaklists;
        this.doNotUseGlobalFDR=doNotUseGlobalFDR;
	}
	
	public void savePreferences(File libraryFile, File fastaFile) throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("encyclopedia");
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
		Preferences prefs=Preferences.userRoot().node("encyclopedia");
		HashMap<String, String> map=new HashMap<String, String>();
		for (String key : prefs.keys()) {
			String value=prefs.get(key, "");
			//System.out.println("Reading EncyclopeDIA preference "+key+" = "+value);
			map.put(key, value);
		}
		return map;
	}

	public String toString() {
		final StringBuilder sb=new StringBuilder();
		sb.append(" -fixed "+aaConstants.getFixedModString()+"\n");
		sb.append(" -frag "+FragmentationType.toString(fragType)+"\n");
		sb.append(" -ptol "+precursorTolerance.getToleranceThreshold()+"\n");
		sb.append(" -ftol "+fragmentTolerance.getToleranceThreshold()+"\n");
		sb.append(" -lftol "+libraryFragmentTolerance.getToleranceThreshold()+"\n");
		sb.append(" -ptolunits "+precursorTolerance.getUnits()+"\n");
		sb.append(" -ftolunits "+fragmentTolerance.getUnits()+"\n");
		sb.append(" -lftolunits "+libraryFragmentTolerance.getUnits()+"\n");
		sb.append(" -poffset "+precursorOffsetPPM+"\n");
		sb.append(" -foffset "+fragmentOffsetPPM+"\n");
		sb.append(" -enzyme "+enzyme.getName()+"\n");
		sb.append(" -percolatorThreshold "+percolatorThreshold+"\n");
		sb.append(" -percolatorVersionNumber "+percolatorVersionNumber+"\n");
		sb.append(" ").append(OPT_PERC_TRAINING_SIZE).append(" ").append(percolatorTrainingSetSize).append("\n");
		sb.append(" ").append(OPT_PERC_TRAINING_THRESH).append(" ").append(percolatorTrainingSetThreshold).append("\n");
		sb.append(" -acquisition "+DataAcquisitionType.toString(dataAcquisitionType)+"\n");
		sb.append(" -numberOfThreadsUsed "+numberOfThreadsUsed+"\n");
		sb.append(" -expectedPeakWidth "+expectedPeakWidth+"\n");
		sb.append(" -precursorWindowSize "+precursorWindowSize+"\n");
		sb.append(" -numberOfQuantitativePeaks "+numberOfQuantitativePeaks+"\n");
		sb.append(" -minNumOfQuantitativePeaks "+minNumOfQuantitativePeaks+"\n");
		sb.append(" -topNTargetsUsed "+topNTargetsUsed+"\n");
		sb.append(" -quantifyAcrossSamples "+quantifyAcrossSamples+"\n");
		sb.append(" -numberOfExtraDecoyLibrariesSearched "+numberOfExtraDecoyLibrariesSearched+"\n");
		sb.append(" -verifyModificationIons "+verifyModificationIons+"\n");
		sb.append(" -minIntensity "+minIntensity+"\n");
		if (useTargetWindowCenter()) {
			sb.append(" -targetWindowCenter "+targetWindowCenter+"\n");
		}
		sb.append(" -scoringBreadthType "+getScoringBreadthType().toShortname()+"\n");
		if (localizingModification.isPresent()) {
			sb.append(" -localizationModification "+localizingModification.get().getShortname()+"\n");
		} else {
			sb.append(" -localizationModification "+PeptideModification.NO_MODIFICATION_NAME+"\n");
		}
		sb.append(" -rtWindowInMin "+rtWindowInMin+"\n");
        sb.append(" -filterPeaklists "+filterPeaklists+"\n");
		return sb.toString();
	}
	
	public HashMap<String, String> toParameterMap() {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("-fixed", aaConstants.getFixedModString());
		map.put("-frag", FragmentationType.toString(fragType));
		map.put("-ptol", precursorTolerance.getToleranceThreshold()+"");
		map.put("-ftol", fragmentTolerance.getToleranceThreshold()+"");
		map.put("-lftol", libraryFragmentTolerance.getToleranceThreshold()+"");
		map.put("-ptolunits", precursorTolerance.getUnits()+"");
		map.put("-ftolunits", fragmentTolerance.getUnits()+"");
		map.put("-lftolunits", libraryFragmentTolerance.getUnits()+"");
		map.put("-poffset", precursorOffsetPPM+"");
		map.put("-foffset", fragmentOffsetPPM+"");
		map.put("-enzyme", enzyme.getName());
		map.put("-percolatorThreshold", percolatorThreshold+"");
		map.put("-percolatorVersionNumber", percolatorVersionNumber+"");
		map.put(OPT_PERC_TRAINING_SIZE, Integer.toString(percolatorTrainingSetSize));
		map.put(OPT_PERC_TRAINING_THRESH, Float.toString(percolatorTrainingSetThreshold));
		map.put("-acquisition", DataAcquisitionType.toString(dataAcquisitionType));
		map.put("-numberOfThreadsUsed", numberOfThreadsUsed+"");
		map.put("-expectedPeakWidth", expectedPeakWidth+"");
		map.put("-precursorWindowSize", precursorWindowSize+"");
		map.put("-numberOfQuantitativePeaks", numberOfQuantitativePeaks+"");
		map.put("-minNumOfQuantitativePeaks", minNumOfQuantitativePeaks+"");
		map.put("-topNTargetsUsed",topNTargetsUsed+"");
		map.put("-quantifyAcrossSamples", quantifyAcrossSamples+"");
		map.put("-numberOfExtraDecoyLibrariesSearched", numberOfExtraDecoyLibrariesSearched+"");
		map.put("-targetWindowCenter", targetWindowCenter+"");
		map.put("-scoringBreadthType", getScoringBreadthType().toShortname());
		map.put("-verifyModificationIons", verifyModificationIons+"");
		map.put("-minIntensity", minIntensity+"");
		if (localizingModification.isPresent()) {
			map.put("-localizationModification", localizingModification.get().getShortname());
		} else {
			map.put("-localizationModification", PeptideModification.NO_MODIFICATION_NAME);
		}
		map.put("-rtWindowInMin", rtWindowInMin+"");
        map.put("-filterPeaklists", filterPeaklists+"");
		return map;
	}
	
	public DataAcquisitionType getDataAcquisitionType() {
		return dataAcquisitionType;
	}
	
	public AminoAcidConstants getAAConstants() {
		return aaConstants;
	}

	public FragmentationType getFragType() {
		return fragType;
	}

	public MassTolerance getFragmentTolerance() {
		return fragmentTolerance;
	}
	
	public MassTolerance getLibraryFragmentTolerance() {
		return libraryFragmentTolerance;
	}

	public MassTolerance getPrecursorTolerance() {
		return precursorTolerance;
	}

	public double getFragmentOffsetPPM() {
		return fragmentOffsetPPM;
	}

	public double getPrecursorOffsetPPM() {
		return precursorOffsetPPM;
	}

	public DigestionEnzyme getEnzyme() {
		return enzyme;
	}

	public float getPercolatorThreshold() {
		return percolatorThreshold;
	}

	public float getPercolatorProteinThreshold() {
		return percolatorProteinThreshold;
	}

	public float getEffectivePercolatorThreshold() {
		return getEffectivePercolatorThreshold(percolatorThreshold, numberOfExtraDecoyLibrariesSearched);
	}

	public static float getEffectivePercolatorThreshold(float percolatorThreshold, float numberOfExtraDecoyLibrariesSearched) {
		// FDR'=FDR * (XD*(1-((XD-1)*FDR)))
		// where XD is the numberOfDecoyLibrariesSearched
		// e.g. if XD=1, then FDR'=FDR*(1*(1-((1-1)*FDR)))=FDR*(1*(1-0))=FDR
		// e.g. if XD=2, then FDR'=FDR*(2*(1-((2-1)*FDR)))=FDR*(2*(1-FDR))=2*FDR-2*FDR*FDR
		float numberOfDecoyLibrariesSearched=numberOfExtraDecoyLibrariesSearched+1.0f; // always search 1x decoy minimum
		return percolatorThreshold*(numberOfDecoyLibrariesSearched*(1-((numberOfDecoyLibrariesSearched-1)*percolatorThreshold)));
	}
	
	public int getPercolatorVersionNumber() {
		return percolatorVersionNumber;
	}

	/**
	 * The value used for the {@code -t/--testFDR} parameter when running
	 * Percolator, which is used to evaluate the best cross-validation
	 * result.
	 *
	 * Currently always {@link PercolatorExecutor#PERCOLATOR_DEFAULT_TEST_THRESHOLD}
	 * to match behavior in version 0.9.4 and earlier. Note that this
	 * differs from the thresholds returned by {@link #getPercolatorThreshold()}
	 * and {@link #getEffectivePercolatorThreshold()} which are applied
	 * to Percolator's results after it's been run. This setting will only
	 * affect the evaluation of cross-validation results, and should not
	 * change the FDR of accepted results.
	 */
	public float getPercolatorTestThreshold() {
		return PercolatorExecutor.PERCOLATOR_DEFAULT_TEST_THRESHOLD;
	}

	/**
	 * The value used for the {@code -N} parameter when running Percolator,
	 * which sets the number of PSMs to use as the training set.
	 */
	public int getPercolatorTrainingSetSize() {
		return percolatorTrainingSetSize;
	}

	/**
	 * The value used for the {@code -F/--trainFDR} parameter when running
	 * Percolator, which sets the FDR used to select positive examples
	 * from the training set.
	 */
	public float getPercolatorTrainingSetThreshold() {
		return percolatorTrainingSetThreshold;
	}

	public boolean isDeconvoluteOverlappingWindows() {
		return dataAcquisitionType==DataAcquisitionType.OVERLAPPING_DIA;
	}
	
	public int getNumberOfThreadsUsed() {
		return numberOfThreadsUsed;
	}
	public float getTargetWindowCenter() {
		return targetWindowCenter;
	}
	public boolean useTargetWindowCenter() {
		return targetWindowCenter>0;
	}
	
	public float getExpectedPeakWidth() {
		return expectedPeakWidth;
	}
	
	public float getPrecursorWindowSize() {
		return precursorWindowSize;
	}
	
	public float getNumberOfExtraDecoyLibrariesSearched() {
		return numberOfExtraDecoyLibrariesSearched;
	}
	
	public int getNumberOfQuantitativePeaks() {
		return numberOfQuantitativePeaks;
	}
	
	public int getEffectiveNumberOfQuantitativePeaks() {
		if (isQuantifySameFragmentsAcrossSamples()) {
			return numberOfQuantitativePeaks;
		} else {
			return Integer.MAX_VALUE;
		}
	}
	
	public int getTopNTargetsUsed() {
		return topNTargetsUsed;
	}
	
	public int getMinNumOfQuantitativePeaks() {
		return minNumOfQuantitativePeaks;
	}
	public float getMinIntensity() {
		return minIntensity;
	}
	public boolean isUseNLsForXCorr() {
		return useNLsForXCorr;
	}
	public double getPrecursorIsolationMargin() {
		return precursorIsolationMargin;
	}
	public ScoringBreadthType getScoringBreadthType() {
		return CASiLBreadthType;
	}
	public Optional<PeptideModification> getLocalizingModification() {
		return localizingModification;
	}
	public boolean applyRTAlignment() {
		return ScoringBreadthType.ENTIRE_RT_WINDOW==CASiLBreadthType;
	}
	
	public boolean isQuantifySameFragmentsAcrossSamples() {
		return quantifyAcrossSamples;
	}
	public boolean isVerifyModificationIons() {
		return verifyModificationIons;
	}
	public float getRtWindowInMin() {
		return rtWindowInMin;
}
    public boolean isFilterPeaklists() {
        return filterPeaklists;
    }
    
    public boolean isDoNotUseGlobalFDR() {
		return doNotUseGlobalFDR;
	}
}