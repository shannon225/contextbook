package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;

public class PercolatorPeptide implements PeptidePrecursor {
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
	
	@Override
	public int hashCode() {
		return getPeptideModSeq().hashCode()+16807*getPrecursorCharge();
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof PeptidePrecursor) {
			return compareTo((PeptidePrecursor)obj)==0;
		}
		return false;
	}
	
	@Override
	public int compareTo(PeptidePrecursor o) {
		if (o==null) return 1;
		int c=getPeptideModSeq().compareTo(o.getPeptideModSeq());
		if (c!=0) return c;
		return Byte.compare(getPrecursorCharge(), o.getPrecursorCharge());
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

	public boolean isPSMIDDecoy() {
		return isPSMIDDecoy(psmID);
	}

	public String getPeptideModSeq() {
		return getPeptideSequence(psmID);
	}

	public byte getPrecursorCharge() {
		return getCharge(psmID);
	}

	public String getFile() {
		return getFile(psmID);
	}

	public float getRT() {
		return getRT(psmID);
	}

	public static String getPSMID(LibraryEntry peptide, float rt, File diaFile) {
		return getPSMID(diaFile.getName(),rt,peptide.isDecoy(),peptide.getPeptideModSeq(),peptide.getPrecursorCharge());
	}
	
	public static String getPSMID(String diaFileName, float rt, boolean isDecoy, String peptideModSeq, byte peptideCharge) {
		return diaFileName+":"+rt+":"+(isDecoy ? "decoy" : "")+peptideModSeq+"+"+peptideCharge;
	}
	

	public static boolean isPSMIDDecoy(String psmID) {
		psmID=getPeptideData(psmID);
		return psmID.startsWith("decoy");
	}

	public static String getPeptideSequence(String psmID) {
		psmID=getPeptideData(psmID);
		if (psmID.startsWith("decoy")) {
			psmID=psmID.substring(5);
		}
		return psmID.substring(0, psmID.lastIndexOf('+'));
	}

	public static byte getCharge(String psmID) {
		return Byte.parseByte(psmID.substring(psmID.lastIndexOf('+')+1));
	}

	public static String getFile(String psmID) {
		int colonIndex=psmID.indexOf(":");
		if (colonIndex>=0) {
			return psmID.substring(0, colonIndex);
		}
		return "";
	}

	public static float getRT(String psmID) {
		int colonIndex=psmID.indexOf(":");
		int colonIndex2=psmID.lastIndexOf(":");
		if (colonIndex2>colonIndex) {
			return Float.parseFloat(psmID.substring(colonIndex+1, colonIndex2));
		}
		return 0.0f;
	}

	public static String getPeptideData(String psmID) {
		int colonIndex=psmID.lastIndexOf(":");
		if (colonIndex>=0) {
			return psmID.substring(colonIndex+1);
		}
		return psmID;
	}
}
