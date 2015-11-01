package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Swath;

public interface PSMScorer {

	public abstract float score(LibraryEntry entry, Swath spectrum);

}