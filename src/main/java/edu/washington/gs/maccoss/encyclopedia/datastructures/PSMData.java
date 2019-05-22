package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.StringTokenizer;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;

public class PSMData implements PeptidePrecursor {
	private final int spectrumIndex;
	private final double precursorMZ;
	private final byte precursorCharge;
	private final String peptideModSeq;
	private final String massCorrectedPeptideModSeq;
	private final float retentionTime;
	private final float score;
	private final float sortingScore;
	private final float duration;
	private final HashSet<String> accessions;
	private final boolean inferred;

	public PSMData(HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, float retentionTime, float score, float sortingScore, float duration, boolean inferred, AminoAcidConstants aaConstants) {
		this(accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, PeptideUtils.getCorrectedMasses(peptideModSeq, aaConstants), retentionTime, score, sortingScore, duration, inferred);
	}
	
	private PSMData(HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, String massCorrectedPeptideModSeq, float retentionTime, float score, float sortingScore, float duration, boolean inferred) {
		this.accessions=accessions;
		this.spectrumIndex=spectrumIndex;
		this.precursorMZ=precursorMZ;
		this.precursorCharge=precursorCharge;
		this.peptideModSeq=peptideModSeq;
		this.massCorrectedPeptideModSeq=massCorrectedPeptideModSeq;
		this.retentionTime=retentionTime;
		this.score=score;
		this.sortingScore=sortingScore;
		this.duration=duration;
		this.inferred=inferred;
	}
	
	public PSMData updateRetentionTime(float rtInSec) {
		return new PSMData(accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, massCorrectedPeptideModSeq, rtInSec, score, sortingScore, duration, inferred);
	}
	
	public boolean wasInferred() {
		return inferred;
	}

	@Override
	public String getPeptideModSeq() {
		return massCorrectedPeptideModSeq;
	}
	
	public double getAlteratelyChargedMass(byte charge) {
		if (charge==getPrecursorCharge()) return getPrecursorMZ();
		
		return (getUnchargedMass()+MassConstants.protonMass*charge)/charge;
	}
	
	public double getUnchargedMass() {
		return getPrecursorMZ()*getPrecursorCharge()-MassConstants.protonMass*getPrecursorCharge();
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

	public String getAccession() {
		ArrayList<String> list=new ArrayList<String>(accessions);
		Collections.sort(list);
		StringBuilder sb=new StringBuilder();
		for (String string : list) {
			if (sb.length()>0) sb.append(PSMData.ACCESSION_TOKEN);
			sb.append(string);
		}
		return sb.toString();
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

	public String getLegacyPeptideModSeq() {
		return peptideModSeq;
	}
	
	public String getPeptideSeq() {
		return PeptideUtils.getPeptideSeq(peptideModSeq);
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
	public static String accessionsToString(Collection<String> accessions) {
		return accessionsToString(accessions, ACCESSION_TOKEN);
	}
	public static String accessionsToString(Collection<String> accessions, String deliminator) {
		if (accessions.size()==0) return "";
		if (accessions.size()==1) return accessions.iterator().next();

		StringBuilder sb=new StringBuilder();
		ArrayList<String> accessionList=new ArrayList<String>(accessions);
		Collections.sort(accessionList);
		for (String string : accessionList) {
			if (sb.length()>0) sb.append(ACCESSION_TOKEN);
			sb.append(string);
		}
		return sb.toString();
	}
	
	public static HashSet<String> stringToAccessions(String string) {
		HashSet<String> accessions=new HashSet<String>();
		if (string!=null) {
			StringTokenizer st = new StringTokenizer(string, ACCESSION_TOKEN);
			while (st.hasMoreTokens()) {
				accessions.add(st.nextToken());
			}
		}
		if (accessions.size()==0) {
			accessions.add("unknown_protein");
		}
		return accessions;
	}
}
