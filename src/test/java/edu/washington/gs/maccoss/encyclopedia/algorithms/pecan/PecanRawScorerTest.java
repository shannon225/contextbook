package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import junit.framework.TestCase;

public class PecanRawScorerTest extends TestCase {

	public void testGetIndividualPeakScores() {
		PecanRawScorer scorer=new PecanRawScorer(new MassTolerance(10.0f), null);
		LibraryEntry entry=getEntry(new double[] {300.0, 1000.01, 1200.0}, new float[] {7, 10, 4});
		Stripe spectrum=getStripe(new double[] {1000.0, 1000.011, 1000.02, 1000.03, 1000.04, 1000.05}, new float[] {10.0f, 100.0f, 50.0f, 10.0f, 1.0f, 100f});
		PeakScores[] scores=scorer.getIndividualPeakScores(entry, spectrum, false);
		
		assertEquals(entry.getMassArray().length, scores.length);
		assertNull(scores[0]);
		assertEquals(scores[1].getScore(), 1600.0, 0.1);
		assertEquals(scores[1].getDeltaMass(), -1, 0.001);
		assertNull(scores[2]);
	}
	
	public LibraryEntry getEntry(double[] masses, float[] intensities) {
		return new LibraryEntry(new HashSet<String>(), 1, (byte)1, "", 1, 1, 1, masses, intensities);
	}
	public Stripe getStripe(double[] masses, float[] intensities) {
		return new Stripe("", "", 1, 1, 1, 1, masses, intensities);
	}
}
