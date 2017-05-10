package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

import gnu.trove.map.hash.TIntIntHashMap;

public class FastaPeptideEntry implements FastaEntryInterface {
	private final String filename;
	private final HashSet<String> accessions=new HashSet<String>();
	private final String sequence;

	public FastaPeptideEntry(String filename, HashSet<String> accessions, String sequence) {
		this.filename=filename;
		this.accessions.addAll(accessions);
		this.sequence=sequence;
	}

	public FastaPeptideEntry(String filename, String accession, String sequence) {
		this.filename=filename;
		accessions.add(accession);
		this.sequence=sequence;
	}

	public FastaPeptideEntry(String sequence) {
		this("Unknown File", "Unknown Annotation", sequence);
	}
	
	public void addAccession(String accession) {
		accessions.add(accession);
	}
	
	public void addAccessions(HashSet<String> accessions) {
		this.accessions.addAll(accessions);
	}

	@Override
	public int compareTo(FastaEntryInterface o) {
		if (o==null) return 1;
		
		int c=getAccession().compareTo(o.getAccession());
		if (c!=0) return c;
		c=filename.compareTo(o.getFilename());
		if (c!=0) return c;
		c=sequence.compareTo(o.getSequence());
		
		return c;
		
	}

	@Override
	public FastaPeptideEntry getSubEntry(String subSequence) {
		return new FastaPeptideEntry(filename, getAccession(), subSequence);
	}
	
	public FastaPeptideEntry getEntryAsPeptide() {
		return new FastaPeptideEntry(filename, accessions, sequence);
	}
	
	public HashSet<String> getAccessions() {
		return accessions;
	}

	@Override
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
	
	public String getAnnotation() {
		return getAccession();
	}

	@Override
	public String getFilename() {
		return filename;
	}

	@Override
	public String getSequence() {
		return sequence;
	}

	@Override
	public void addStatistics(TIntIntHashMap map) {
		FastaEntry.getStatistics(sequence, map);
	}

}
