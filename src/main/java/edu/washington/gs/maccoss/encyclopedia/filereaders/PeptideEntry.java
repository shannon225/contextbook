package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

class PeptideEntry {
	final String peptideModSeq;
	final float rt;
	final byte charge;
	final ArrayList<Peak> peaks;

	public PeptideEntry(String peptideModSeq, byte charge, float rt) {
		this.peptideModSeq=peptideModSeq;
		this.charge=charge;
		this.rt=rt;
		this.peaks=new ArrayList<>();
	}

	public void addPeak(Peak peak) {
		peaks.add(peak);
	}

	@Override
	public String toString() {
		return peptideModSeq+"+"+charge+","+rt+" ("+peaks.size()+"): "+General.toString(peaks);
	}
}