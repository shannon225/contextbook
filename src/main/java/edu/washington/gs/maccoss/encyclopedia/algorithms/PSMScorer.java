package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;

public interface PSMScorer {

	public abstract float score(LibraryEntry entry, Stripe spectrum);

}