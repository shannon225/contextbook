package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.Quadruplet;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class QuantitativeDIAData implements PeptidePrecursor, Spectrum {

	private final String peptideModSeq;
	private final String massCorrectedPeptideModSeq;
	private final byte precursorCharge;
	private final double precursorMZ;
	private final float scanStartTime;
	private final Range rtScanRange;
	private final double[] massArray;
	private final float[] intensityArray;

	public QuantitativeDIAData(String peptideModSeq, byte precursorCharge, float scanStartTime, Range rtScanRange, double[] massArray, float[] intensityArray, AminoAcidConstants aaConstants) {
		this.peptideModSeq=peptideModSeq;
		this.massCorrectedPeptideModSeq=PeptideUtils.getCorrectedMasses(peptideModSeq, aaConstants);
		this.precursorMZ=aaConstants.getChargedMass(massCorrectedPeptideModSeq, precursorCharge);
		this.precursorCharge=precursorCharge;
		this.scanStartTime=scanStartTime;
		
		ArrayList<PeakChromatogram> peaks=new ArrayList<>();
		int numPeaks=massArray.length;
		for (int i=0; i<numPeaks; i++) {
			peaks.add(new PeakChromatogram(massArray[i], intensityArray[i], 0.0f, false));
		}
		Collections.sort(peaks);
		Quadruplet<double[], float[], float[], boolean[]> arrays=PeakChromatogram.toChromatogramArrays(peaks);
		this.massArray=arrays.x;
		this.intensityArray=arrays.y;
		
		this.rtScanRange=rtScanRange;
	}
	
	@Override
	public float getScanStartTime() {
		return scanStartTime;
	}
	@Override
	public double getPrecursorMZ() {
		return precursorMZ;
	}
	@Override
	public String getSpectrumName() {
		return getPeptideModSeq();
	}

	@Override
	public String getPeptideModSeq() {
		return massCorrectedPeptideModSeq;
	}

	@Override
	public int compareTo(PeptidePrecursor o) {
		if (o==null) return 1;
		int c=getPeptideModSeq().compareTo(o.getPeptideModSeq());
		if (c!=0) return c;
		return Byte.compare(getPrecursorCharge(), o.getPrecursorCharge());
	}
	
	@Override
	public int hashCode() {
		return getPeptideModSeq().hashCode()+16807*getPrecursorCharge();
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof PeptidePrecursor) {
			return compareTo((PeptidePrecursor)obj)==0;
		}
		return false;
	}

	public String getLegacyPeptideModSeq() {
		return peptideModSeq;
	}

	public String getPeptideSeq() {
		StringBuilder sb=new StringBuilder();
		for (char c : getPeptideModSeq().toCharArray()) {
			if (Character.isLetter(c)) {
				sb.append(c);
			}
		}
		return sb.toString();
	}
	
	public Range getRtScanRange() {
		return rtScanRange;
	}

	public byte getPrecursorCharge() {
		return precursorCharge;
	}

	public float getApexRT() {
		return scanStartTime;
	}

	public double[] getMassArray() {
		return massArray;
	}

	public float[] getIntensityArray() {
		return intensityArray;
	}

	public float getTIC() {
		return General.sum(intensityArray);
	}
	
	public int getNumNonZeroPeaks() {
		int n=0;
		for (int i=0; i<intensityArray.length; i++) {
			if (intensityArray[i]>0.0f) n++;
		}
		return n;
	}

}
