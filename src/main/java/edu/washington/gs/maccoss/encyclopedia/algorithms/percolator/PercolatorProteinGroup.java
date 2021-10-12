package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroupInterface;

public class PercolatorProteinGroup implements ProteinGroupInterface {
	private final String[] accessions;
	private final String[] peptides;
	private final float qValue;
	private final float posteriorErrorProb;
	private final int hash;

	public PercolatorProteinGroup(List<String> accessions, List<String> peptides, float qValue, float posteriorErrorProb) {
		this(accessions.toArray(new String[accessions.size()]), peptides.toArray(new String[peptides.size()]), qValue, posteriorErrorProb);
	}
	
	public PercolatorProteinGroup(String[] accessions, String[] peptides, float qValue, float posteriorErrorProb) {
		this.accessions=accessions;
		Arrays.sort(this.accessions);
		
		this.peptides=peptides;
		Arrays.sort(this.peptides);
		
		this.qValue=qValue;
		this.posteriorErrorProb=posteriorErrorProb;

		hash=PSMData.accessionsToString(this.getEquivalentAccessions()).hashCode();
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
		return PSMData.accessionsToString(getEquivalentAccessions());
	}

	@Override
	public int compareTo(ProteinGroupInterface o) {
		if (o==null) return 1;

		int c=Float.compare(getNSPScore(), o.getNSPScore());
		if (c!=0) return c;

		List<String> acc2=o.getEquivalentAccessions();
		List<String> acc1=getEquivalentAccessions();
		c=Integer.compare(acc1.size(), acc2.size());
		if (c!=0) return c;
		return PSMData.accessionsToString(acc1).compareTo(PSMData.accessionsToString(acc2));
	}

	@Override
	public float getNSPScore() {
		// good enough estimation
		return (1.0f-posteriorErrorProb)*peptides.length;
	}

	@Override
	public List<String> getEquivalentAccessions() {
		return Arrays.asList(accessions);
	}

	@Override
	public List<String> getSequences() {
		return Arrays.asList(peptides);
	}

	public HashSet<String> getAccessions() {
		return new HashSet<>(Arrays.asList(accessions));
	}

	public String[] getPeptides() {
		return peptides;
	}

	public float getQValue() {
		return qValue;
	}

	public float getPosteriorErrorProb() {
		return posteriorErrorProb;
	}

	@Override
	public boolean isDecoy() {
		for (String accession : accessions) {
			if (!accession.startsWith(LibraryEntry.DECOY_STRING)) {
				return false;
			}
		}
		return true;
	}

}
