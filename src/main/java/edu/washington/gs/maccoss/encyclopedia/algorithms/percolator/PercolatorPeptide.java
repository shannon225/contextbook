package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;

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

	public boolean isPSMIDDecoy() {
		return isPSMIDDecoy(psmID);
	}

	public String getPeptideSequence() {
		return getPeptideSequence(psmID);
	}

	public byte getCharge() {
		return getCharge(psmID);
	}

	public String getFile() {
		return getFile(psmID);
	}

	public float getRT() {
		return getRT(psmID);
	}

	public static String getPSMID(LibraryEntry peptide, float rt, File diaFile) {
		return diaFile.getName()+":"+rt+":"+(peptide.isDecoy() ? "decoy" : "")+peptide.getPeptideModSeq()+"+"+peptide.getPrecursorCharge();
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
