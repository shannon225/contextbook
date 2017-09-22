package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;
import java.util.Iterator;
import java.util.stream.Collectors;

import com.google.common.collect.ImmutableSortedSet;
import gnu.trove.map.hash.TIntIntHashMap;

public class FastaPeptideEntry implements Comparable<FastaPeptideEntry> {
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
	public int compareTo(FastaPeptideEntry o) {
		if (o==null) return 1;
		
		int c = getAccessions().size() - o.getAccessions().size();
		if (c!=0) return c;

		final Iterator<String>
				iter = ImmutableSortedSet.copyOf(getAccessions()).iterator(),
				oIter = ImmutableSortedSet.copyOf(o.getAccessions()).iterator();
		while (iter.hasNext()) {
			String acc = iter.next(), oAcc = oIter.next();

			c = acc.compareTo(oAcc);
			if (c!=0) return c;
		}

		c=filename.compareTo(o.getFilename());
		if (c!=0) return c;
		c=sequence.compareTo(o.getSequence());
		
		return c;
		
	}

	public HashSet<String> getAccessions() {
		return accessions;
	}

	/**
	 * Returns this entry's set of accessions with the given flag
	 * prepended to each. Useful for generating decoys.
	 */
	public HashSet<String> getFlaggedAccessions(String prefix) {
		return getAccessions().stream()
				.map(s -> prefix + s)
				.collect(Collectors.toCollection(HashSet::new));
	}

	public String getFilename() {
		return filename;
	}

	public String getSequence() {
		return sequence;
	}

	public void addStatistics(TIntIntHashMap map) {
		FastaEntry.getStatistics(sequence, map);
	}

}
