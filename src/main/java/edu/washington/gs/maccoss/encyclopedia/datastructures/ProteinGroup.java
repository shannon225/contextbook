package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Collections;

public class ProteinGroup implements Comparable<ProteinGroup> {
	private final float nspScore;
	private final ArrayList<String> equivalentAccessions;
	private final int hash;

	/**
	 * 
	 * @param nspScore
	 * @param equivalentAccessions note, destructively sorts this array!
	 */
	public ProteinGroup(float nspScore, ArrayList<String> equivalentAccessions) {
		this.nspScore=nspScore;
		this.equivalentAccessions=equivalentAccessions;
		Collections.sort(equivalentAccessions);
		
		hash=getAccessionString(equivalentAccessions).hashCode();
	}

	private String getAccessionString(ArrayList<String> equivalentAccessions) {
		return PSMData.accessionsToString(equivalentAccessions);
	}
	
	@Override
	public int hashCode() {
		return hash;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj==null||!(obj instanceof ProteinGroup)) return false;
		if (hashCode()!=obj.hashCode()) return false;
		return equivalentAccessions.toString().equals(((ProteinGroup)obj).equivalentAccessions.toString());
	}
	
	@Override
	public String toString() {
		return getAccessionString(equivalentAccessions);
	}
	
	@Override
	public int compareTo(ProteinGroup o) {
		if (o==null) return 1;
		int c=Float.compare(nspScore, o.nspScore);
		if (c!=0) return c;
		
		c=Integer.compare(equivalentAccessions.size(), o.equivalentAccessions.size());
		if (c!=0) return c;
		return equivalentAccessions.toString().compareTo(((ProteinGroup)o).equivalentAccessions.toString());
	}

	public float getNspScore() {
		return nspScore;
	}

	public ArrayList<String> getEquivalentAccessions() {
		return equivalentAccessions;
	}
}
