package edu.washington.gs.maccoss.encyclopedia.filereaders;

import gnu.trove.map.hash.TIntIntHashMap;

public class FastaEntry {
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
	
	public void addStatistics(TIntIntHashMap map) {
		getStatistics(sequence, map);
	}

	public static void getStatistics(String sequence, TIntIntHashMap map) {
		for (int i=0; i<sequence.length(); i++) {
			map.adjustOrPutValue(sequence.charAt(i), 1, 1);
		}
	}
}
