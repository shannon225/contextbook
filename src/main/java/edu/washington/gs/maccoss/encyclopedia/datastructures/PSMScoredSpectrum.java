package edu.washington.gs.maccoss.encyclopedia.datastructures;

public class PSMScoredSpectrum implements HasMZ, HasRetentionTime, HasCharge {
	private final PSMData peptideData;
	private final FragmentScan msms;

	public PSMScoredSpectrum(PSMData peptideData, FragmentScan msms) {
		this.peptideData = peptideData;
		this.msms = msms;
	}

	public FragmentScan getMsms() {
		return msms;
	}

	public PSMData getPeptideData() {
		return peptideData;
	}

	@Override
	public float getRetentionTimeInSec() {
		return peptideData.getRetentionTimeInSec();
	}

	@Override
	public double getMZ() {
		return peptideData.getPrecursorMZ();
	}
	
	@Override
	public byte getCharge() {
		return peptideData.getPrecursorCharge();
	}
	
}
