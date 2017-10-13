package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import junit.framework.TestCase;

public class EncyclopediaOneScorerTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false);

	public void testGetIndividualPeakScores() {
		
		EncyclopediaOneScorer scorer=new EncyclopediaOneScorer(PARAMETERS, null);
		LibraryEntry entry=getEntry(new double[] {300.0, 1000.01, 1200.0}, new float[] {7, 10, 4});
		Stripe spectrum=getStripe(new double[] {1000.0, 1000.011, 1000.02, 1000.03, 1000.04, 1000.05}, new float[] {10.0f, 100.0f, 50.0f, 10.0f, 1.0f, 100f});
		FragmentIon[] targets=new FragmentIon[] {
				new FragmentIon(300, (byte)2, IonType.y),
				new FragmentIon(1000.01, (byte)4, IonType.y),
				new FragmentIon(1200, (byte)5, IonType.y)
		};
		
		PeakScores[] scores=scorer.getIndividualPeakScores(entry, spectrum, false, targets);
		for (int i=0; i<scores.length; i++) {
			if (scores[i]!=null) {
				System.out.println(scores[i].getTargetMass()+"\t"+scores[i].getScore()+"\t"+scores[i].getDeltaMass());
			} else {
				System.out.println("null");
			}
		}
		
		assertEquals(targets.length, scores.length);
		assertNull(scores[0]);
		assertEquals(scores[1].getScore(), 1600.0, 0.1);
		assertEquals(scores[1].getDeltaMass(), -1, 0.001);
		assertNull(scores[2]);
	}
	
	public LibraryEntry getEntry(double[] masses, float[] intensities) {
		return new LibraryEntry("", new HashSet<String>(), 1, (byte)1, "", 1, 1, 1, masses, intensities);
	}
	public Stripe getStripe(double[] masses, float[] intensities) {
		return new Stripe("", "", 1, 1, 1, 1, masses, intensities);
	}
}
