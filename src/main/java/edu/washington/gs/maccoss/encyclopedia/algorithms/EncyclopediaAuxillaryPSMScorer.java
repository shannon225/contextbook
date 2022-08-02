package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public abstract class EncyclopediaAuxillaryPSMScorer extends AuxillaryPSMScorer {

	public EncyclopediaAuxillaryPSMScorer(SearchParameters parameters) {
		super(parameters);
	}

	public EncyclopediaAuxillaryPSMScorer getEntryOptimizedScorer(LibraryEntry entry) {
		throw new EncyclopediaException("getEntryOptimizedScorer must be overridden!");
	}
}
