package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import edu.washington.gs.maccoss.encyclopedia.Encyclopedia;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public class SearchParameters {

	protected final AminoAcidConstants aaConstants;
	protected final FragmentationType fragType;
	protected final MassTolerance precursorTolerance;
	protected final MassTolerance fragmentTolerance;
	protected final MassTolerance libraryFragmentTolerance;
	protected final DigestionEnzyme enzyme;
	protected final float percolatorThreshold;
	protected final File percolatorLocation;
	protected final DataAcquisitionType dataAcquisitionType;
	protected final int numberOfThreadsUsed;	
	protected final float targetWindowCenter;
	protected final float expectedPeakWidth;
	protected final boolean runPhosphoLocalization;
	protected final float precursorWindowSize;
	protected final float numberOfExtraDecoyLibrariesSearched;
	protected final int numberOfQuantitativePeaks;
	protected final double precursorOffsetPPM;
	protected final double fragmentOffsetPPM;

	public SearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance precursorTolerance, double precursorOffsetPPM, MassTolerance fragmentTolerance, double fragmentOffsetPPM, MassTolerance libraryFragmentTolerance, DigestionEnzyme enzyme,
			float percolatorThreshold, File percolatorLocation, DataAcquisitionType dataAcquisitionType, int numberOfThreadsUsed, float expectedPeakWidth, float targetWindowCenter, float precursorWindowSize, 
			int numberOfQuantitativePeaks, boolean runPhosphoLocalization, float getNumberOfExtraDecoyLibrariesSearched) {
		this.aaConstants=aaConstants;
		this.fragType=fragType;
		this.precursorTolerance=precursorTolerance;
		this.precursorOffsetPPM=precursorOffsetPPM;
		this.fragmentTolerance=fragmentTolerance;
		this.fragmentOffsetPPM=fragmentOffsetPPM;
		this.libraryFragmentTolerance=libraryFragmentTolerance;
		this.enzyme=enzyme;
		this.percolatorThreshold=percolatorThreshold;
		this.percolatorLocation=percolatorLocation;
		this.dataAcquisitionType=dataAcquisitionType;
		this.numberOfThreadsUsed=numberOfThreadsUsed;
		this.expectedPeakWidth=expectedPeakWidth;
		this.targetWindowCenter=targetWindowCenter;
		this.precursorWindowSize=precursorWindowSize;
		this.numberOfQuantitativePeaks=numberOfQuantitativePeaks;
		this.runPhosphoLocalization=runPhosphoLocalization;
		this.numberOfExtraDecoyLibrariesSearched=getNumberOfExtraDecoyLibrariesSearched;
	}
	
	public void savePreferences(File libraryFile) throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("encyclopedia");
		HashMap<String, String> map=toParameterMap();
		if (libraryFile!=null) map.put(Encyclopedia.TARGET_LIBRARY_TAG, libraryFile.getAbsolutePath());
		for (Entry<String, String> entry : map.entrySet()) {
			//System.out.println("Writing EncyclopeDIA preference "+entry.getKey()+" = "+entry.getValue());
			prefs.put(entry.getKey(), entry.getValue());
		}
		prefs.flush();
		//OutputStream stream=new FileOutputStream("encyclopedia.prefs");
		//prefs.exportNode(stream);
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
		sb.append(" -ptol "+precursorTolerance.getPpmTolerance()+"\n");
		sb.append(" -ftol "+fragmentTolerance.getPpmTolerance()+"\n");
		sb.append(" -lftol "+libraryFragmentTolerance.getPpmTolerance()+"\n");
		sb.append(" -poffset "+precursorOffsetPPM+"\n");
		sb.append(" -foffset "+fragmentOffsetPPM+"\n");
		sb.append(" -enzyme "+enzyme.getName()+"\n");
		sb.append(" -percolatorThreshold "+percolatorThreshold+"\n");
		sb.append(" -percolatorLocation "+percolatorLocation+"\n");
		sb.append(" -acquisition "+DataAcquisitionType.toString(dataAcquisitionType)+"\n");
		sb.append(" -numberOfThreadsUsed "+numberOfThreadsUsed+"\n");
		sb.append(" -expectedPeakWidth "+expectedPeakWidth+"\n");
		sb.append(" -precursorWindowSize "+precursorWindowSize+"\n");
		sb.append(" -numberOfQuantitativePeaks "+numberOfQuantitativePeaks+"\n");
		sb.append(" -runPhosphoLocalization "+runPhosphoLocalization+"\n");
		sb.append(" -getNumberOfExtraDecoyLibrariesSearched "+numberOfExtraDecoyLibrariesSearched+"\n");
		if (useTargetWindowCenter()) {
			sb.append(" -targetWindowCenter "+targetWindowCenter+"\n");
		}
		return sb.toString();
	}
	
	public HashMap<String, String> toParameterMap() {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("-fixed", aaConstants.getFixedModString());
		map.put("-frag", FragmentationType.toString(fragType));
		map.put("-ptol", precursorTolerance.getPpmTolerance()+"");
		map.put("-ftol", fragmentTolerance.getPpmTolerance()+"");
		map.put("-lftol", libraryFragmentTolerance.getPpmTolerance()+"");
		map.put("-poffset", precursorOffsetPPM+"");
		map.put("-foffset", fragmentOffsetPPM+"");
		map.put("-enzyme", enzyme.getName());
		map.put("-percolatorThreshold", percolatorThreshold+"");
		map.put("-percolatorLocation", percolatorLocation+"");
		map.put("-acquisition", DataAcquisitionType.toString(dataAcquisitionType));
		map.put("-numberOfThreadsUsed", numberOfThreadsUsed+"");
		map.put("-expectedPeakWidth", expectedPeakWidth+"");
		map.put("-precursorWindowSize", precursorWindowSize+"");
		map.put("-numberOfQuantitativePeaks", numberOfQuantitativePeaks+"");
		map.put("-runPhosphoLocalization", runPhosphoLocalization+"");
		map.put("-getNumberOfExtraDecoyLibrariesSearched", numberOfExtraDecoyLibrariesSearched+"");
		map.put("-targetWindowCenter", targetWindowCenter+"");
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

	public float getEffectivePercolatorThreshold() {
		// FDR'=FDR * (XD*(1-((XD-1)*FDR)))
		// where XD is the numberOfDecoyLibrariesSearched
		// e.g. if XD=1, then FDR'=FDR*(1*(1-((1-1)*FDR)))=FDR*(1*(1-0))=FDR
		// e.g. if XD=2, then FDR'=FDR*(2*(1-((2-1)*FDR)))=FDR*(2*(1-FDR))=2*FDR-2*FDR*FDR
		float numberOfDecoyLibrariesSearched=numberOfExtraDecoyLibrariesSearched+1.0f; // always search 1x decoy minimum
		return percolatorThreshold*(numberOfDecoyLibrariesSearched*(1-((numberOfDecoyLibrariesSearched-1)*percolatorThreshold)));
	}

	public Optional<File> getPercolatorLocation() {
		return Optional.ofNullable(percolatorLocation);
	}

	public boolean isDeconvoluteOverlappingWindows() {
		return dataAcquisitionType==DataAcquisitionType.OVERLAPPING_DIA;
	}
	
	public boolean isDDA() {
		return dataAcquisitionType==DataAcquisitionType.DDA;
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
	
	public boolean isRunPhosphoLocalization() {
		return runPhosphoLocalization;
	}
	public float getNumberOfExtraDecoyLibrariesSearched() {
		return numberOfExtraDecoyLibrariesSearched;
	}
	public int getNumberOfQuantitativePeaks() {
		return numberOfQuantitativePeaks;
	}
}