package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ProteinGroup implements ProteinGroupInterface {
	private final float nspScore;
	private final float posteriorErrorProb;
	private final ArrayList<String> equivalentAccessions;
	private final int hash;
	private final ArrayList<String> sequences;

	/**
	 * 
	 * @param nspScore
	 * @param equivalentAccessions note, destructively sorts this array!
	 */
	public ProteinGroup(float nspScore, float posteriorErrorProb, ArrayList<String> equivalentAccessions, ArrayList<String> sequences) {
		this.nspScore=nspScore;
		this.posteriorErrorProb=posteriorErrorProb;
		this.equivalentAccessions=equivalentAccessions;
		this.sequences=sequences;
		Collections.sort(this.equivalentAccessions);
		
		hash=PSMData.accessionsToString(this.equivalentAccessions).hashCode();
	}
	
	@Override
	public int hashCode() {
		return hash;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj==null||!(obj instanceof ProteinGroupInterface)) return false;
		if (hashCode()!=obj.hashCode()) return false;
		return PSMData.accessionsToString(getEquivalentAccessions()).equals(PSMData.accessionsToString(((ProteinGroupInterface)obj).getEquivalentAccessions()));
	}
	
	@Override
	public String toString() {
		return PSMData.accessionsToString(equivalentAccessions);
	}
	
	@Override
	public ArrayList<String> getSequences() {
		return sequences;
	}

	@Override
	public float getNSPScore() {
		return nspScore;
	}
	
	@Override
	public int compareTo(ProteinGroupInterface o) {
		if (o==null) return 1;
		
		//int c=Float.compare(getNSPScore(), o.getNSPScore());
		int c=-Float.compare(getPosteriorErrorProb(), o.getPosteriorErrorProb());
		if (c!=0) return c;
		
		List<String> acc2=o.getEquivalentAccessions();
		List<String> acc1=getEquivalentAccessions();
		c=Integer.compare(acc1.size(), acc2.size());
		if (c!=0) return c;
		return PSMData.accessionsToString(acc1).compareTo(PSMData.accessionsToString(acc2));
	}

	@Override
	public ArrayList<String> getEquivalentAccessions() {
		return equivalentAccessions;
	}
	
	@Override
	public boolean isDecoy() {
		for (String accession : equivalentAccessions) {
			if (!accession.startsWith(LibraryEntry.DECOY_STRING)) {
				return false;
			}
		}
		return true;
	}
	
	@Override
	public float getPosteriorErrorProb() {
		return posteriorErrorProb;
	}
}
