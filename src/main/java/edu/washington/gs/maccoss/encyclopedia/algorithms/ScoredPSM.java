package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeptideXYPoint;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;

public class ScoredPSM implements Comparable<ScoredPSM>, PSMInterface {
	private final LibraryEntry entry;
	private final FragmentScan msms;
	private final float primaryScore;
	private final float[] auxScores;
	private final float deltaPrecursorMass;
	private final float deltaFragmentMass;
	
	public ScoredPSM(LibraryEntry entry, FragmentScan msms, float primaryScore, float[] auxScores,
			float deltaPrecursorMass, float deltaFragmentMass) {
		this.entry = entry;
		this.msms = msms;
		this.primaryScore = primaryScore;
		this.auxScores = auxScores;
		this.deltaPrecursorMass = deltaPrecursorMass;
		this.deltaFragmentMass = deltaFragmentMass;
	}

	public LibraryEntry getLibraryEntry() {
		return entry;
	}

	public FragmentScan getMSMS() {
		return msms;
	}

	public float getPrimaryScore() {
		return primaryScore;
	}

	public float[] getAuxScores() {
		return auxScores;
	}
	
	public float getDeltaFragmentMass() {
		return deltaFragmentMass;
	}
	
	public float getDeltaPrecursorMass() {
		return deltaPrecursorMass;
	}
	
	public PeptideXYPoint getRTData() {
		return new PeptideXYPoint(entry.getScanStartTime()/60f, msms.getScanStartTime()/60f, entry.isDecoy(), entry.getPeptideModSeq());
	}

	@Override
	public int compareTo(ScoredPSM o) {
		if (o==null) return 1;
		int c=Float.compare(primaryScore, o.primaryScore);
		if (c!=0) return c;
		c=entry.compareTo(o.entry);
		if (c!=0) return c;
		c=msms.compareTo(o.msms);
		if (c!=0) return c;
		return 0;
	}

	@Override
	public int hashCode() {
		return Float.hashCode(primaryScore);
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof ScoredPSM) return compareTo((ScoredPSM)obj)==0;
		return false;
	}
}
