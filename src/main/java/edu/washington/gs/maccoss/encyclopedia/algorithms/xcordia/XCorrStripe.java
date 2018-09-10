package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrSpectrum;

public class XCorrStripe extends Stripe {
	private final SparseXCorrSpectrum xcorrSpectrum;
	
	public XCorrStripe(Stripe stripe, SearchParameters params) {
		super(stripe.getSpectrumName(), stripe.getPrecursorName(), stripe.getSpectrumIndex(), stripe.getScanStartTime(), stripe.getIonInjectionTime(), stripe.getIsolationWindowLower(), stripe.getIsolationWindowUpper(), stripe.getMassArray(), stripe.getIntensityArray());
		xcorrSpectrum=SparseXCorrCalculator.normalize(this, new Range(stripe.getIsolationWindowLower(), stripe.getIsolationWindowUpper()), false, params);
	}

	public SparseXCorrSpectrum getXcorrSpectrum() {
		return xcorrSpectrum;
	}
	
	public static ArrayList<Stripe> downcastXCorrToStripe(ArrayList<XCorrStripe> stripes) {
		ArrayList<Stripe> downcast=new ArrayList<Stripe>();
		for (Stripe stripe : stripes) {
			downcast.add(stripe);
		}
		return downcast;
	}
}
