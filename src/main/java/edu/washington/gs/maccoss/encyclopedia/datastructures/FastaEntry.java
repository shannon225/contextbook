package edu.washington.gs.maccoss.encyclopedia.datastructures;

import gnu.trove.map.hash.TIntIntHashMap;

public class FastaEntry implements FastaEntryInterface {
	private final String filename;
	private final String annotation;
	private final String sequence;

	public FastaEntry(String sequence) {
		this("Unknown File", "Unknown Annotation", sequence);
	}
	public FastaEntry(String filename, String annotation, String sequence) {
		this.filename=filename;
		this.annotation=annotation;
		this.sequence=sequence;
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntryInterface#getAccession()
	 */
	@Override
	public String getAccession() {
		int start=annotation.charAt(0)=='>'?1:0;
		int stop=annotation.indexOf(' ');
		if (stop<0) {
			stop=annotation.length();
		}
		return annotation.substring(start, stop);
	}
	
	public String getAnnotation() {
		int start=annotation.charAt(0)=='>'?1:0;
		return annotation.substring(start);
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntryInterface#getFilename()
	 */
	@Override
	public String getFilename() {
		return filename;
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntryInterface#getSequence()
	 */
	@Override
	public String getSequence() {
		return sequence;
	}

	@Override
	public FastaPeptideEntry getSubEntry(String subSequence) {
		return new FastaPeptideEntry(filename, getAccession(), subSequence);
	}
	
	public FastaPeptideEntry getEntryAsPeptide() {
		return new FastaPeptideEntry(filename, getAccession(), sequence);
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntryInterface#addStatistics(gnu.trove.map.hash.TIntIntHashMap)
	 */
	@Override
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
	public int compareTo(FastaEntryInterface o) {
		if (o==null) return 1;
		
		int c=getAccession().compareTo(o.getAccession());
		if (c!=0) return c;
		c=filename.compareTo(o.getFilename());
		if (c!=0) return c;
		c=sequence.compareTo(o.getSequence());
		
		return c;
		
	}
	
}
