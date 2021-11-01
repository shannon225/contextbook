package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.QuantitativeDIAData;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;

public class PeptideReportData implements PeptidePrecursor {
	private final String massCorrectedPeptideModSeq;
	private final byte precursorCharge;
	private final String accessions;
	
	volatile private float avgRT=-1;
	volatile private int maxNumOfFragments=0;

	public PeptideReportData(String peptideModSeq, byte precursorCharge, String accessions, AminoAcidConstants aaConstants) {
		this.massCorrectedPeptideModSeq=PeptideUtils.getCorrectedMasses(peptideModSeq, aaConstants);
		this.precursorCharge=precursorCharge;
		this.accessions=accessions;
	}
	
	HashMap<String, QuantitativeDIAData> dataBySourceFile=new HashMap<>();
	public void addQuantitativeDIAData(String sourceFile, QuantitativeDIAData data) {
		final QuantitativeDIAData previousData = dataBySourceFile.put(sourceFile, data);
		if (null != previousData) {
			throw new IllegalStateException("Data conflicts with existing entry for sample " + sourceFile + " (charge " + data.getPrecursorCharge() + " vs. " + previousData.getPrecursorCharge() + ")");
		}
		maxNumOfFragments=Math.max(maxNumOfFragments, data.getMassArray().length);
	}
	
	public int getMaxNumOfFragments() {
		return maxNumOfFragments;
	}

	public byte getPrecursorCharge() {
		return precursorCharge;
	}

	public String getAccessionString() {
		return accessions;
	}
	
	public HashSet<String> getAccessions() {
		return PSMData.stringToAccessions(accessions);
	}

	public QuantitativeDIAData getQuantitativeData(String sourceFile) {
		return dataBySourceFile.get(sourceFile);
	}
	
	public Set<String> getSampleNames() {
		return dataBySourceFile.keySet();
	}

	public float getAverageRetentionTime() {
		if (avgRT<0) {
			TFloatArrayList rtCenters=new TFloatArrayList();
			for (QuantitativeDIAData data : dataBySourceFile.values()) {
				rtCenters.add(data.getApexRT());
			}
			avgRT=General.mean(rtCenters.toArray());
		}
		return avgRT;
	}

	@Override
	public String getPeptideModSeq() {
		return massCorrectedPeptideModSeq;
	}

	@Override
	public int compareTo(PeptidePrecursor o) {
		if (o==null) return 1;
		int c=getPeptideModSeq().compareTo(o.getPeptideModSeq());
		if (c!=0) return c;
		return Byte.compare(getPrecursorCharge(), o.getPrecursorCharge());
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

//	@Override
//	public String getLegacyPeptideModSeq() {
//		return peptideModSeq;
//	}

	@Override
	public String getPeptideSeq() {
		StringBuilder sb=new StringBuilder();
		for (char c : getPeptideModSeq().toCharArray()) {
			if (Character.isLetter(c)) {
				sb.append(c);
			}
		}
		return sb.toString();
	}
}
