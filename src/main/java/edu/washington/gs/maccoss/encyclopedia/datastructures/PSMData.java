package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;
import java.util.StringTokenizer;

public class PSMData implements PeptidePrecursor {
	private final int spectrumIndex;
	private final double precursorMZ;
	private final byte precursorCharge;
	private final String peptideModSeq;
	private final float retentionTime;
	private final float score;
	private final float sortingScore;
	private final float duration;
	private final HashSet<String> accessions;

	public PSMData(HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, float retentionTime, float score, float sortingScore, float duration) {
		this.accessions=accessions;
		this.spectrumIndex=spectrumIndex;
		this.precursorMZ=precursorMZ;
		this.precursorCharge=precursorCharge;
		this.peptideModSeq=peptideModSeq;
		this.retentionTime=retentionTime;
		this.score=score;
		this.sortingScore=sortingScore;
		this.duration=duration;
	}
	
	@Override
	public int compareTo(PeptidePrecursor o) {
		if (o==null) return 1;
		int c=getPeptideModSeq().compareTo(o.getPeptideModSeq());
		if (c!=0) return c;
		return Byte.compare(getPrecursorCharge(), o.getPrecursorCharge());
	}
	
	public HashSet<String> getAccessions() {
		return accessions;
	}

	public int getSpectrumIndex() {
		return spectrumIndex;
	}

	public double getPrecursorMZ() {
		return precursorMZ;
	}

	public byte getPrecursorCharge() {
		return precursorCharge;
	}

	public String getPeptideModSeq() {
		return peptideModSeq;
	}

	public float getRetentionTime() {
		return retentionTime;
	}

	public float getScore() {
		return score;
	}
	
	public float getSortingScore() {
		return sortingScore;
	}
	
	public float getDuration() {
		return duration;
	}
	
	public static final String ACCESSION_TOKEN=";";
	public static String accessionsToString(HashSet<String> accessions) {
		StringBuilder sb=new StringBuilder();
		for (String string : accessions) {
			if (sb.length()>0) sb.append(ACCESSION_TOKEN);
			sb.append(string);
		}
		return sb.toString();
	}
	
	public static HashSet<String> stringToAccessions(String string) {
		StringTokenizer st=new StringTokenizer(string, ACCESSION_TOKEN);
		HashSet<String> accessions=new HashSet<String>();
		while (st.hasMoreTokens()) {
			accessions.add(st.nextToken());
		}
		if (accessions.size()==0) {
			accessions.add("unknown_protein");
		}
		return accessions;
	}
}
