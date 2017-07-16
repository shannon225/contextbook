package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.Collection;

public class PeptideAccessionMatchingTrie extends PeptideTrie<LibraryEntry> {
	
	public PeptideAccessionMatchingTrie(Collection<LibraryEntry> entries) {
		super(entries);
	}

	@Override
	protected void processMatch(FastaEntryInterface fasta, LibraryEntry entry, int start) {
		entry.getAccessions().add(fasta.getAccession());
	}
}
