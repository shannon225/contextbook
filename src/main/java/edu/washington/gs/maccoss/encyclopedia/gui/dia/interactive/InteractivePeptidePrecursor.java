package edu.washington.gs.maccoss.encyclopedia.gui.dia.interactive;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.HasRetentionTime;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SimplePeptidePrecursor;

public class InteractivePeptidePrecursor extends SimplePeptidePrecursor implements HasRetentionTime {
	private static final AminoAcidConstants aaConstants=new AminoAcidConstants();
	private final float rtInSecs;
	private final double precursorMZ;
	
	private byte isPassing; // CAN BE MUTATED!
	private Range rtRangeInSecs; // CAN BE MUTATED!

	public InteractivePeptidePrecursor(String peptideModSeq, byte precursorCharge, float rtInSecs) {
		this(peptideModSeq, precursorCharge, rtInSecs, (byte)0);
	}

	public InteractivePeptidePrecursor(String peptideModSeq, byte precursorCharge, float rtInSecs, byte isPassing) {
		super(peptideModSeq, precursorCharge, aaConstants);
		this.rtInSecs=rtInSecs;
		this.precursorMZ=aaConstants.getChargedMass(peptideModSeq, precursorCharge);
		this.isPassing=isPassing;
	}

	@Override
	public float getRetentionTimeInSec() {
		return rtInSecs;
	}
	
	public double getPrecursorMZ() {
		return precursorMZ;
	}
	
	/**
	 * -1 for not passing
	 *  0 for unassigned
	 *  1 for passing
	 * @return
	 */
	public byte getIsPassing() {
		return isPassing;
	}

	/**
	 * -1 for not passing
	 *  0 for unassigned
	 *  1 for passing
	 */
	public void setIsPassing(boolean pass) {
		if (pass) {
			isPassing=1;
		} else {
			isPassing=-1;
		}
	}
	
	public Range getRTRange() {
		return rtRangeInSecs;
	}
	public void setRtRangeInSecs(Range rtRangeInSecs) {
		this.rtRangeInSecs=rtRangeInSecs;
	}
}
