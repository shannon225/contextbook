package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

public class ScoredMProphetData implements Comparable<ScoredMProphetData> {
	MProphetData data;
	float score;
	double pvalue;
	double localFDR;
	double fdr;
	public ScoredMProphetData(MProphetData data, float score, double pvalue, double localFDR, double fdr) {
		this.data = data;
		this.score = score;
		this.pvalue = pvalue;
		this.localFDR = localFDR;
		this.fdr = fdr;
	}
	
	public MProphetData getData() {
		return data;
	}
	public double getFDR() {
		return fdr;
	}
	public double getLocalFDR() {
		return localFDR;
	}
	public double getPvalue() {
		return pvalue;
	}
	public float getScore() {
		return score;
	}
	
	@Override
	public int compareTo(ScoredMProphetData o) {
		if (o==null) return -1;
		int c=Double.compare(localFDR, o.localFDR);
		if (c!=0) return c;
		
		return data.compareTo(o.data);
	}
}
