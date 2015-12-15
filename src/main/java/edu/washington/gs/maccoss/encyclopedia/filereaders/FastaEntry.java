package edu.washington.gs.maccoss.encyclopedia.filereaders;

import gnu.trove.map.hash.TIntIntHashMap;

public class FastaEntry implements Comparable<FastaEntry> {
	private final String filename;
	private final String annotation;
	private final String sequence;

	public FastaEntry(String filename, String annotation, String sequence) {
		this.filename=filename;
		this.annotation=annotation;
		this.sequence=sequence;
	}

	public String getAnnotation() {
		return annotation;
	}

	public String getFilename() {
		return filename;
	}

	public String getSequence() {
		return sequence;
	}
	
	public FastaEntry getSubEntry(String subSequence) {
		return new FastaEntry(filename, annotation, subSequence);
	}
	
	public void addStatistics(TIntIntHashMap map) {
		getStatistics(sequence, map);
	}

	public static void getStatistics(String sequence, TIntIntHashMap map) {
		for (int i=0; i<sequence.length(); i++) {
			map.adjustOrPutValue(sequence.charAt(i), 1, 1);
		}
	}
	
	@Override
	public int hashCode() {
		int hashCode=1;
		hashCode=31*hashCode+(filename==null?0:filename.hashCode());
		hashCode=31*hashCode+(annotation==null?0:annotation.hashCode());
		hashCode=31*hashCode+(sequence==null?0:sequence.hashCode());
		return hashCode;
	}
	
	@Override
	public boolean equals(Object obj) {
		return compareTo((FastaEntry)obj)==0;
	}
	@Override
	public int compareTo(FastaEntry o) {
		if (o==null) return 1;
		
		int c=annotation.compareTo(o.annotation);
		if (c!=0) return c;
		c=filename.compareTo(o.filename);
		if (c!=0) return c;
		c=sequence.compareTo(o.sequence);
		
		return c;
		
	}
	
}
