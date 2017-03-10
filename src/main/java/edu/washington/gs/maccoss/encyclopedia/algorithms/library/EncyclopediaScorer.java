package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;

public interface EncyclopediaScorer extends PSMPeakScorer {

	EncyclopediaOneAuxillaryPSMScorer getAuxScorer();

	float score(LibraryEntry entry, Spectrum spectrum);

	float score(LibraryEntry entry, Spectrum spectrum, FragmentIon[] ions);

}