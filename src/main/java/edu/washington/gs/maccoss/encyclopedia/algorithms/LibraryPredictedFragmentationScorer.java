package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import gnu.trove.list.array.TFloatArrayList;

public class LibraryPredictedFragmentationScorer extends AuxillaryPSMScorer {

	public LibraryPredictedFragmentationScorer(SearchParameters parameters) {
		super(parameters);
	}

	@Override
	public float[] score(LibraryEntry entry, Stripe spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		MassTolerance tolerance=parameters.getFragmentTolerance();
		double[] masses=spectrum.getMassArray();
		float[] intensities=spectrum.getIntensityArray();
		double[] predicted=entry.getMassArray();
		
		TFloatArrayList ions=new TFloatArrayList();
		for (int i=0; i<predicted.length; i++) {
			int[] indicies=tolerance.getIndicies(masses, predicted[i]);
			float intensity=0.0f;
			for (int j=0; j<indicies.length; j++) {
				intensity+=intensities[indicies[j]];
			}
			ions.add(intensity);
		}
		return ions.toArray();
	}

	@Override
	public float[] getMissingDataScores(LibraryEntry entry) {
		return new float[entry.getMassArray().length];
	}

	@Override
	public String[] getScoreNames(LibraryEntry entry) {
		ArrayList<String> names=new ArrayList<String>();

		double[] predicted=entry.getMassArray();
		for (int i=0; i<predicted.length; i++) {
			names.add(Long.toString(Math.round(predicted[i])));
		}
		return names.toArray(new String[names.size()]);
	}
}
