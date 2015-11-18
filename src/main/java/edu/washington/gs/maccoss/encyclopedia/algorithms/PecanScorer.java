package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import jdk.nashorn.internal.ir.annotations.Immutable;

@Immutable
public class PecanScorer implements PSMScorer {
	private final MassTolerance fragmentTolerance;
	private final MassTolerance precursorTolerance;
	private final PrecursorScanMap precursors;


	public PecanScorer(MassTolerance fragmentTolerance, MassTolerance precursorTolerance, PrecursorScanMap precursors) {
		this.fragmentTolerance=fragmentTolerance;
		this.precursorTolerance=precursorTolerance;
		this.precursors=precursors;
	}


	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer#score(edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry, edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe)
	 */
	public float score(LibraryEntry entry, Stripe spectrum) {
		// precursor scoring
		Peak[] precursorPacket=precursors.getIsotopePacket(entry.getPrecursorMZ(), spectrum.getScanStartTime(), entry.getPrecursorCharge(), precursorTolerance);
		
		// FIXME
		// calculate mass accuracy
		// calculate isotope dot product (with -1 peak)
		
		
		
		// fragment scoring
		double[] libraryMasses=entry.getMassArray();
		float[] libraryIntensities=entry.getIntensityArray();
		
		double[] spectrumMasses=spectrum.getMassArray();
		float[] spectrumIntensities=spectrum.getIntensityArray();
		
		if (libraryMasses.length==0||spectrumMasses.length==0) return 0.0f;
		
		float sum=0.0f;
		int libraryIndex=0;
		int spectrumIndex=0;
		while (true) {
			int compare=fragmentTolerance.compareTo(libraryMasses[libraryIndex], spectrumMasses[spectrumIndex]);
			if (compare==0) {
				sum+=libraryIntensities[libraryIndex]*spectrumIntensities[spectrumIndex];
				libraryIndex++;
				spectrumIndex++;
			} else if (compare>0) {
				spectrumIndex++;
			} else {
				libraryIndex++;
			}
			if (libraryIndex>=libraryMasses.length) break;
			if (spectrumIndex>=spectrumMasses.length) break;
		}
		
		return sum;
	}
}
