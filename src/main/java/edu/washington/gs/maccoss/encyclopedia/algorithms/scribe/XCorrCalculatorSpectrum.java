package edu.washington.gs.maccoss.encyclopedia.algorithms.scribe;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorrStripe;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrSpectrum;

public class XCorrCalculatorSpectrum extends FragmentScan {
	private final SparseXCorrCalculator xcorrCalculator;
	private final int numPeaks;
	
	public XCorrCalculatorSpectrum(FragmentScan stripe, SearchParameters params) {
		super(stripe.getSpectrumName(), stripe.getPrecursorName(), stripe.getSpectrumIndex(), stripe.getScanStartTime(), stripe.getFraction(), stripe.getIonInjectionTime(), stripe.getIsolationWindowLower(), stripe.getIsolationWindowUpper(), stripe.getMassArray(), stripe.getIntensityArray());

		SparseXCorrSpectrum xcorrSpectrum=SparseXCorrCalculator.normalize(stripe, new Range(stripe.getIsolationWindowLower(), stripe.getIsolationWindowUpper()), false, params);
		numPeaks=xcorrSpectrum.getIndices().length;
		this.xcorrCalculator=new SparseXCorrCalculator(xcorrSpectrum, params);
		
	}
	
	public int getNumPeaks() {
		return numPeaks;
	}
	
	public static ArrayList<FragmentScan> downcastXCorrToStripe(ArrayList<XCorrStripe> stripes) {
		ArrayList<FragmentScan> downcast=new ArrayList<FragmentScan>();
		for (FragmentScan stripe : stripes) {
			downcast.add(stripe);
		}
		return downcast;
	}

	public float score(SparseXCorrSpectrum entry) {
		return xcorrCalculator.score(entry);
	}
}
