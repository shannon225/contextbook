package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.util.Comparator;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;

public class PercolatorPeptide implements PeptidePrecursor {
	private final String psmID;
	private final String proteinIDs;
	private final float qValue;
	private final float posteriorErrorProb;
	private final String massCorrectedPeptideModSeq;
	
	public static final Comparator<PercolatorPeptide> scoreComparator=new Comparator<PercolatorPeptide>() {
		// lower is better
		@Override
		public int compare(PercolatorPeptide o1, PercolatorPeptide o2) {
			if (o1==null&&o2==null) return 0;
			if (o1==null) return -1;
			if (o2==null) return 1;
			
			int c=Float.compare(o1.qValue, o2.qValue);
			if (c!=0) return c;
			c=Float.compare(o1.posteriorErrorProb, o2.posteriorErrorProb);
			if (c!=0) return c;
			return o1.compareTo(o2);
		}
	};

	public PercolatorPeptide(String psmID, String proteinIDs, float qValue, float posteriorErrorProb) {
		this.psmID=psmID;
		this.proteinIDs=proteinIDs;
		this.qValue=qValue;
		this.posteriorErrorProb=posteriorErrorProb;
		this.massCorrectedPeptideModSeq=PeptideUtils.getCorrectedMasses(getPeptideSequence(psmID));
	}
	
	@Override
	public String getMassCorrectedPeptideModSeq() {
		return massCorrectedPeptideModSeq;
	}
	
	public String getPeptideSeq() {
		StringBuilder sb=new StringBuilder();
		for (char c : getPeptideModSeq().toCharArray()) {
			if (Character.isLetter(c)) {
				sb.append(c);
			}
		}
		return sb.toString();
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
		int c=getMassCorrectedPeptideModSeq().compareTo(o.getMassCorrectedPeptideModSeq());
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

	public static String getPSMID(LibraryEntry peptide, float rt, StripeFileInterface diaFile) {
		return getPSMID(diaFile.getOriginalFileName(),rt,peptide.isDecoy(),peptide.getPeptideModSeq(),peptide.getPrecursorCharge());
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
