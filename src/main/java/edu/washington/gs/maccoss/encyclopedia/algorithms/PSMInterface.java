package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeptideXYPoint;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;

public interface PSMInterface {
	public PeptidePrecursor getLibraryEntry();

	public float getDeltaFragmentMass();
	
	public float getDeltaPrecursorMass();
	
	public PeptideXYPoint getRTData();
}
