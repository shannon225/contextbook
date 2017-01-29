package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrSpectrum;

public class XCorrStripe extends Stripe {
	private final SparseXCorrSpectrum xcorrSpectrum;


	public XCorrStripe(String spectrumName, String precursorName, int spectrumIndex, float scanStartTime, float isolationWindowLower, float isolationWindowUpper, double[] massArray,
			float[] intensityArray, SearchParameters params) {
		super(spectrumName, precursorName, spectrumIndex, scanStartTime, isolationWindowLower, isolationWindowUpper, massArray, intensityArray);
		xcorrSpectrum=SparseXCorrCalculator.normalize(this, new Range(isolationWindowLower, isolationWindowUpper), false, params);
	}

	public SparseXCorrSpectrum getXcorrSpectrum() {
		return xcorrSpectrum;
	}
	
	public static ArrayList<Stripe> downcast(ArrayList<XCorrStripe> stripes) {
		ArrayList<Stripe> downcast=new ArrayList<Stripe>();
		for (Stripe stripe : stripes) {
			downcast.add(stripe);
		}
		return downcast;
	}
}
