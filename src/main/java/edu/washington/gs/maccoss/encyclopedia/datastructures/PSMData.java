package edu.washington.gs.maccoss.encyclopedia.datastructures;

public class PSMData {
	private final int spectrumIndex;
	private final double precursorMZ;
	private final byte precursorCharge;
	private final String peptideModSeq;
	private final float retentionTime;
	private final float score;
	private final float duration;

	public PSMData(int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, float retentionTime, float score, float duration) {
		this.spectrumIndex=spectrumIndex;
		this.precursorMZ=precursorMZ;
		this.precursorCharge=precursorCharge;
		this.peptideModSeq=peptideModSeq;
		this.retentionTime=retentionTime;
		this.score=score;
		this.duration=duration;
	}

	public int getSpectrumIndex() {
		return spectrumIndex;
	}

	public double getPrecursorMZ() {
		return precursorMZ;
	}

	public byte getPrecursorCharge() {
		return precursorCharge;
	}

	public String getPeptideModSeq() {
		return peptideModSeq;
	}

	public float getRetentionTime() {
		return retentionTime;
	}

	public float getScore() {
		return score;
	}
	
	public float getDuration() {
		return duration;
	}
}
