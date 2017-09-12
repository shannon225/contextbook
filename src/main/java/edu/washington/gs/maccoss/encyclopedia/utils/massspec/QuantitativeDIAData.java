package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class QuantitativeDIAData implements PeptidePrecursor {

	private final String peptideModSeq;
	private final String massCorrectedPeptideModSeq;
	private final byte precursorCharge;
	private final float scanStartTime;
	private final double[] massArray;
	private final float[] intensityArray;

	public QuantitativeDIAData(String peptideModSeq, byte precursorCharge, float scanStartTime, double[] massArray, float[] intensityArray) {
		this.peptideModSeq=peptideModSeq;
		this.massCorrectedPeptideModSeq=PeptideUtils.getCorrectedMasses(peptideModSeq);
		this.precursorCharge=precursorCharge;
		this.scanStartTime=scanStartTime;
		this.massArray=massArray;
		this.intensityArray=intensityArray;
	}

	@Override
	public String getMassCorrectedPeptideModSeq() {
		return massCorrectedPeptideModSeq;
	}

	@Override
	public int compareTo(PeptidePrecursor o) {
		if (o==null) return 1;
		int c=getMassCorrectedPeptideModSeq().compareTo(o.getMassCorrectedPeptideModSeq());
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

	public String getPeptideModSeq() {
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

	public byte getPrecursorCharge() {
		return precursorCharge;
	}

	public float getScanStartTime() {
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
