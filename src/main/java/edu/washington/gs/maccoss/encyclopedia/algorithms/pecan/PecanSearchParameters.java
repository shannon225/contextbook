package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;

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
		sb.append(" -ptol "+precursorTolerance.getPpmTolerance()+"\n");
		sb.append(" -ftol "+fragmentTolerance.getPpmTolerance()+"\n");
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
		sb.append(" -deconvoluteOverlappingWindows "+dataAcquisitionType+"\n");
		sb.append(" -numberOfThreadsUsed "+numberOfThreadsUsed+"\n");
		sb.append(" -alpha "+alpha+"\n");
		sb.append(" -beta "+beta+"\n");
		
		return sb.toString();
	}
	
	public PecanSearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance precursorTolerance, MassTolerance fragmentTolerance, DigestionEnzyme enzyme, int minPeptideLength,
			int maxPeptideLength, int maxMissedCleavages, byte minCharge, byte maxCharge, int minEluteTime, int numberOfReportedPeaks, boolean addDecoysToBackgound, boolean dontRunDecoys, float percolatorThreshold, float alpha, float beta, File percolatorLocation, DataAcquisitionType dataAcquisitionType, int numberOfThreadsUsed, float targetWindowCenter, float precursorWindowSize) {
		super(aaConstants, fragType, precursorTolerance, fragmentTolerance, enzyme, percolatorThreshold, percolatorLocation, dataAcquisitionType, numberOfThreadsUsed, minEluteTime*2.0f, targetWindowCenter, precursorWindowSize, false);
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
	
	public PecanSearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance precursorTolerance, MassTolerance fragmentTolerance, DigestionEnzyme enzyme,
			int maxMissedCleavages, byte minCharge, byte maxCharge, DataAcquisitionType dataAcquisitionType, float precursorWindowSize, int numberOfJobs) {
		super(aaConstants, fragType, precursorTolerance, fragmentTolerance, enzyme, 0.01f, null, dataAcquisitionType, numberOfJobs, 24f, -1f, precursorWindowSize, false);
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

	public PecanSearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance fragmentTolerance, MassTolerance precursorTolerance, DigestionEnzyme enzyme) {
		super(aaConstants, fragType, precursorTolerance, fragmentTolerance, enzyme, 0.01f, null, DataAcquisitionType.DIA, Runtime.getRuntime().availableProcessors(), 24f, -1f, -1f, false);
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
			int maxMissedCleavages) {
		super(aaConstants, fragType, precursorTolerance, fragmentTolerance, enzyme, 0.01f, null, DataAcquisitionType.DIA, Runtime.getRuntime().availableProcessors(), 24f, -1f, -1f, false);
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
