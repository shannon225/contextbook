package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class QuantitativeDIAData implements PeptidePrecursor {

	private final String peptideModSeq;
	private final String massCorrectedPeptideModSeq;
	private final byte precursorCharge;
	private final float scanStartTime;
	private final Range rtScanRange;
	private final double[] massArray;
	private final float[] intensityArray;

	public QuantitativeDIAData(String peptideModSeq, byte precursorCharge, float scanStartTime, Range rtScanRange, double[] massArray, float[] intensityArray, AminoAcidConstants aaConstants) {
		this.peptideModSeq=peptideModSeq;
		this.massCorrectedPeptideModSeq=PeptideUtils.getCorrectedMasses(peptideModSeq, aaConstants);
		this.precursorCharge=precursorCharge;
		this.scanStartTime=scanStartTime;
		this.massArray=massArray;
		this.intensityArray=intensityArray;
		this.rtScanRange=rtScanRange;
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
