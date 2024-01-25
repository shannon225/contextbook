package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.util.Collections;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;

public class ImmutablePeptideEntry implements Comparable<ImmutablePeptideEntry> {
	final String sourceFile;
	final String peptideModSeq;
	final float rt;
	final byte charge;
	final double[] masses;
	final float[] intensities;
	
	public ImmutablePeptideEntry(PeptideEntry entry) {
		peptideModSeq=entry.peptideModSeq;
		rt=entry.rt;
		charge=entry.charge;
		sourceFile=entry.sourceFile;
		
		Collections.sort(entry.peaks);
		Triplet<double[], float[], Optional<float[]>> peakArrays=Peak.toArrays(entry.peaks);
		masses=peakArrays.x;
		intensities=peakArrays.y;
	}
	
	@Override
	public int compareTo(ImmutablePeptideEntry o) {
		if (o==null) return 1;
		int c=this.charge-o.charge;
		if (c!=0) return c;
		c=this.peptideModSeq.compareTo(o.peptideModSeq);
		return c;
	}
}
