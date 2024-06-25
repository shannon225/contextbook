package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeptideXYPoint;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;

public class SimplePSM implements PSMInterface {
	private final PeptidePrecursor entry;
	private final float deltaPrecursorMass;
	private final float deltaFragmentMass;
	private final PeptideXYPoint rtPoint;
	
	public SimplePSM(PeptidePrecursor entry, float deltaPrecursorMass, float deltaFragmentMass,
			PeptideXYPoint rtPoint) {
		this.entry = entry;
		this.deltaPrecursorMass = deltaPrecursorMass;
		this.deltaFragmentMass = deltaFragmentMass;
		this.rtPoint = rtPoint;
	}

	@Override
	public PeptidePrecursor getLibraryEntry() {
		// TODO Auto-generated method stub
		return entry;
	}

	@Override
	public float getDeltaFragmentMass() {
		return deltaFragmentMass;
	}

	@Override
	public float getDeltaPrecursorMass() {
		return deltaPrecursorMass;
	}

	@Override
	public PeptideXYPoint getRTData() {
		return rtPoint;
	}

}
