package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

public class PercolatorPeptide {
	private final String psmID;
	private final String proteinIDs;
	private final float qValue;
	private final float posteriorErrorProb;

	public PercolatorPeptide(String psmID, String proteinIDs, float qValue, float posteriorErrorProb) {
		this.psmID=psmID;
		this.proteinIDs=proteinIDs;
		this.qValue=qValue;
		this.posteriorErrorProb=posteriorErrorProb;
	}

	public String getPsmID() {
		return psmID;
	}

	public String getProteinIDs() {
		return proteinIDs;
	}

	public float getQValue() {
		return qValue;
	}

	public float getPosteriorErrorProb() {
		return posteriorErrorProb;
	}

}
