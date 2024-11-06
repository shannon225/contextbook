package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import junit.framework.TestCase;
import org.junit.Assert;

import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: caleb
 * Date: 10/31/2024
 * Time: 8:53 AM
 */
public class SparseXCorrCalculatorTest extends TestCase {
	// all that matters about these parameters is their tolerances
	private static final SearchParameters MAIN_PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID,
			new MassTolerance(1.0, MassErrorUnitType.AMU), new MassTolerance(1.0, MassErrorUnitType.AMU), DigestionEnzyme.getEnzyme("trypsin"), false, true, false);


	public void testNormalize() {
		Spectrum s = XCorrCalculatorTest.getSpectrum(new double[] { 400.0, 409.0, 421.0, 437.0, 451.0, 498.0, 500.0},
				//                 									410, 410,  430, 440, 460 (X), MAX, MAX // the upper bin
														new float[] { 10f, 2f, 25f, 1f,  1f,      10f, 5f },
				 57f,
				46.7f,
				"toy spectrum",
				451.0
				);

		Range precursorMz = new Range(450.0, 452.0);

		SparseXCorrSpectrum normalized = SparseXCorrCalculator.normalize(s, precursorMz, false, MAIN_PARAMETERS);

		Assert.assertArrayEquals(new double[] {400.0, 409.0, 421.0, 437.0, 498.0, 500.0}, // 451 removed (precursor mz)
				normalized.getMassArray(), 0.00005);
		final float pi = ArrayXCorrCalculator.primaryIonIntensity;
		Assert.assertArrayEquals(new float[] {pi, pi/5, pi, pi, pi, pi/2 },
				normalized.getIntensityArray(), 0.00005f);
	}


	public void testGetTheoreticalSpectrum() {
		final SparseXCorrSpectrum lmk = SparseXCorrCalculator.getTheoreticalSpectrum("LM[+16.0]K", (byte) 2, MAIN_PARAMETERS);
		final float pi = ArrayXCorrCalculator.primaryIonIntensity;

		Assert.assertArrayEquals(new double[] {
				114.0913404668129, 147.1128041505165, 261.12674046681286,
				294.14820415051645, 389.22170346681287, 407.23226815051646}, lmk.getMassArray(), 0.00005);
		Assert.assertArrayEquals(new float[] {pi,pi,pi,pi,pi,pi}, lmk.getIntensityArray(), 0.00005f);
	}
}