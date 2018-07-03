package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import gnu.trove.list.array.TFloatArrayList;

public class ExpectedFragmentationScorer extends AuxillaryPSMScorer {
	private final int startIonIndex; //1=start with b2, y2, etc
	

	public ExpectedFragmentationScorer(SearchParameters parameters, int startIonIndex) {
		super(parameters);
		this.startIonIndex=startIonIndex;
	}

	@Override
	public float[] score(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		MassTolerance tolerance=parameters.getFragmentTolerance();
		double[] masses=spectrum.getMassArray();
		float[] intensities=spectrum.getIntensityArray();
		

		FragmentIon[] modelIons=model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), false);

		TFloatArrayList ions=new TFloatArrayList();
		for (FragmentIon ion : modelIons) {
			if (ion.index<startIonIndex) continue;
			int[] indicies=tolerance.getIndicies(masses, ion.mass);
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
		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		FragmentIon[] modelIons=model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), false);

		TFloatArrayList ions=new TFloatArrayList();
		for (FragmentIon ion : modelIons) {
			if (ion.index<startIonIndex) continue;
			ions.add(0.0f);
		}
		
		return ions.toArray();
	}

	@Override
	public String[] getScoreNames(LibraryEntry entry) {
		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		FragmentIon[] ions=model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), false);
		
		ArrayList<String> names=new ArrayList<String>();
		for (FragmentIon ion : ions) {
			if (ion.index<startIonIndex) continue;
			names.add(ion.toString());
		}
		return names.toArray(new String[names.size()]);
	}
}
