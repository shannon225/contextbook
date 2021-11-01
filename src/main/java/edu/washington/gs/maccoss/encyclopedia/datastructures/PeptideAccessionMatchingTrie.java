package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.Collection;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;

public class PeptideAccessionMatchingTrie extends PeptideTrie<LibraryEntry> {

	public PeptideAccessionMatchingTrie(Collection<LibraryEntry> entries) {
		super(entries);
	}
	
	public PeptideAccessionMatchingTrie(Collection<LibraryEntry> entries, Optional<DigestionEnzyme> enzyme) {
		super(entries, enzyme);
	}

	@Override
	protected void processMatch(FastaEntryInterface fasta, LibraryEntry entry, int start) {
		entry.getAccessions().add(fasta.getAccession());
	}
}
