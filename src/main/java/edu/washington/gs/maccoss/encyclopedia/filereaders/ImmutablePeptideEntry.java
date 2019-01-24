package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;

public class ImmutablePeptideEntry {
	final String peptideModSeq;
	final float rt;
	final byte charge;
	final double[] masses;
	final float[] intensities;
	
	public ImmutablePeptideEntry(PeptideEntry entry) {
		peptideModSeq=entry.peptideModSeq;
		rt=entry.rt;
		charge=entry.charge;
		
		Collections.sort(entry.peaks);
		Pair<double[], float[]> peakArrays=Peak.toArrays(entry.peaks);
		masses=peakArrays.x;
		intensities=peakArrays.y;
	}
}
