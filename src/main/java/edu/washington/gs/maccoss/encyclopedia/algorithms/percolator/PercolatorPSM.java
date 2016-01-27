package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.util.ArrayList;

public class PercolatorPSM implements Comparable<PercolatorPSM> {
	private final String psmID;
	private final boolean isDecoy;
	private final float svmScore;
	private final float qValue;
	private final float pepValue;
	private final float pValue;
	private final String sequence;
	private final char nTermChar;
	private final char cTermChar;
	private final String proteinID;

	public PercolatorPSM(String psmID, boolean isDecoy, float svmScore, float qValue, float pepValue, float pValue, String sequence, char nTermChar, char cTermChar, String proteinID) {
		this.psmID=psmID;
		this.isDecoy=isDecoy;
		this.svmScore=svmScore;
		this.qValue=qValue;
		this.pepValue=pepValue;
		this.pValue=pValue;
		this.sequence=sequence;
		this.nTermChar=nTermChar;
		this.cTermChar=cTermChar;
		this.proteinID=proteinID;
	}
	
	@Override
	public int compareTo(PercolatorPSM o) {
		if (o==null) return 1;
		int c=Float.compare(svmScore, o.svmScore);
		if (c!=0) return c;
		return psmID.compareTo(o.psmID);
	}
	
	public String toPeptideIndex() {
		return isDecoy?"decoy":""+sequence;
	}

	public String toPSMString(float adjustedQValue) {
		return     
		"    <psm p:psm_id=\""+psmID+"\" p:decoy=\""+isDecoy+"\">\n"+
	    "      <svm_score>"+svmScore+"</svm_score>\n"+
	    "      <q_value>"+adjustedQValue+"</q_value>\n"+
	    "      <pep>"+pepValue+"</pep>\n"+
	    "      <exp_mass>"+0.0+"</exp_mass>\n"+
	    "      <calc_mass>"+0.0+"</calc_mass>\n"+
	    "      <peptide_seq n=\""+nTermChar+"\" c=\""+cTermChar+"\" seq=\""+sequence+"\"/>\n"+
	    "      <protein_id>"+proteinID+"</protein_id>\n"+
	    "      <p_value>"+pValue+"</p_value>\n"+
	    "    </psm>\n";
	}
	
	public static String toPeptideString(ArrayList<PercolatorPSM> psms, float adjustedQValue) {
		PercolatorPSM bestPSM=null;
		StringBuilder sb=new StringBuilder();
		for (PercolatorPSM percolatorPSM : psms) {
			sb.append("        <psm_id>"+percolatorPSM.psmID+"</psm_id>\n");
			
			if (bestPSM==null||percolatorPSM.svmScore>bestPSM.svmScore) {
				bestPSM=percolatorPSM;
			} else if (percolatorPSM.svmScore==bestPSM.svmScore) { 
				// don't have to null check since nulls get swallowed into the above logic
				if (!percolatorPSM.isDecoy&&bestPSM.isDecoy) {
					// prefer non-decoys if equal
					bestPSM=percolatorPSM;
				}
			}
		}
		
		if (bestPSM==null) return "";
		
		return     
		"    <peptide p:peptide_id=\""+bestPSM.sequence+"\" p:decoy=\""+bestPSM.isDecoy+"\">\n"+
	    "      <svm_score>"+bestPSM.svmScore+"</svm_score>\n"+
	    "      <q_value>"+adjustedQValue+"</q_value>\n"+
	    "      <pep>"+bestPSM.pepValue+"</pep>\n"+
	    "      <exp_mass>"+0.0+"</exp_mass>\n"+
	    "      <calc_mass>"+0.0+"</calc_mass>\n"+
	    "      <protein_id>"+bestPSM.proteinID+"</protein_id>\n"+
	    "      <p_value>"+bestPSM.pValue+"</p_value>\n"+
	    "      <psm_ids>\n"+
	    sb.toString()+
	    "      </psm_ids>\n"+
	    "    </peptide>\n";
	}
}
