package edu.washington.gs.maccoss.encyclopedia.gui.dia.interactive;

import java.util.StringTokenizer;

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
		this(peptideModSeq, precursorCharge, rtInSecs, isPassing, null);
	}

	public InteractivePeptidePrecursor(String peptideModSeq, byte precursorCharge, float rtInSecs, byte isPassing, Range rtRangeInSecs) {
		super(peptideModSeq, precursorCharge, aaConstants);
		this.rtInSecs=rtInSecs;
		this.precursorMZ=aaConstants.getChargedMass(peptideModSeq, precursorCharge);
		this.isPassing=isPassing;
		this.rtRangeInSecs=rtRangeInSecs;
	}
	
	public static String getHeader() {
		return "#Peptide\tRT(min)\tCharge\tIsPassing\tRTStart(min)\tRTStop(min)";
	}
	
	public String toString() {
		StringBuilder sb=new StringBuilder();
		sb.append(getPeptideModSeq());
		sb.append('\t');
		sb.append(rtInSecs/60f);
		sb.append('\t');
		sb.append(getPrecursorCharge());
		sb.append('\t');
		sb.append(getIsPassing());
		
		if (rtRangeInSecs==null) {
			sb.append('\t');
			sb.append('\t');
		} else {
			sb.append('\t');
			sb.append(rtRangeInSecs.getStart()/60f);
			sb.append('\t');
			sb.append(rtRangeInSecs.getStop()/60f);
		}
		return sb.toString();
	}
	
	public static InteractivePeptidePrecursor parseFromLine(String line) {
		StringTokenizer st = new StringTokenizer(line);
		String peptideModSeq = st.nextToken();
		float rtMin = Float.parseFloat(st.nextToken());
		byte charge = Byte.parseByte(st.nextToken());
		byte passes = st.hasMoreTokens() ? Byte.parseByte(st.nextToken()) : 0;
		
		Range range=null;
		if (st.hasMoreTokens()) {
			float first =Float.parseFloat(st.nextToken());
			float second =Float.parseFloat(st.nextToken());
			range=new Range(first*60f, second*60f);
		}

		return new InteractivePeptidePrecursor(peptideModSeq, charge, rtMin * 60f, passes, range);
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
	public void removeIsPassing() {
		isPassing=0;
	}
	
	/**
	 * CAN BE NULL
	 * @return
	 */
	public Range getRTRange() {
		return rtRangeInSecs;
	}
	public void setRtRangeInSecs(Range rtRangeInSecs) {
		this.rtRangeInSecs=rtRangeInSecs;
	}
}
