package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;

public interface PSMPeakScorer extends PSMScorer {
	public abstract PeakScores[] getIndividualPeakScores(LibraryEntry entry, Spectrum spectrum, boolean normalize);
	public abstract PeakScores[] getIndividualPeakScores(LibraryEntry entry, Spectrum spectrum, boolean normalize, FragmentIon[] ions);
}