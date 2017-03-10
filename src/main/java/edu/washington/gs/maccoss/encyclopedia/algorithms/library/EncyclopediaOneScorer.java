package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;

public class EncyclopediaOneScorer implements EncyclopediaScorer {
	private final SearchParameters parameters;
	private final EncyclopediaOneAuxillaryPSMScorer auxScorer;

	public EncyclopediaOneScorer(SearchParameters parameters, LibraryBackground background) {
		this.parameters=parameters;
		auxScorer=new EncyclopediaOneAuxillaryPSMScorer(parameters, background, true);
	}
	

	@Override
	public EncyclopediaOneAuxillaryPSMScorer getAuxScorer() {
		return auxScorer;
	}
	
	@Override
	public float[] auxScore(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return auxScorer.score(entry, spectrum, predictedIsotopeDistribution, precursors);
	}
	@Override
	public String[] getAuxScoreNames(LibraryEntry entry) {
		return auxScorer.getScoreNames(entry);
	}

	@Override
	public float score(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return score(entry, spectrum);
	}

	@Override
	public float score(LibraryEntry entry, Spectrum spectrum) {
		PeakScores[] individualPeakScores=getIndividualPeakScores(entry, spectrum, true);
		return scoreIons(individualPeakScores);
	}


	@Override
	public float score(LibraryEntry entry, Spectrum spectrum, FragmentIon[] ions) {
		PeakScores[] individualPeakScores=getIndividualPeakScores(entry, spectrum, true, ions);
		return scoreIons(individualPeakScores);
	}

	private float scoreIons(PeakScores[] individualPeakScores) {
		int count=0; // number of matches
		for (int i=0; i<individualPeakScores.length; i++) {
			if (individualPeakScores[i]!=null) count++;
		}
		
		if (count==0) return 0.0f;
		
		float dotProduct=PeakScores.sumScores(individualPeakScores); // dot product
		
		return Log.protectedLog10(dotProduct)+Log.logFactorial(count); // X!Tandem score
	}

	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Spectrum spectrum, boolean normalize) {
		FragmentationModel model=new FragmentationModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		FragmentIon[] ions=model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge());
		
		return getIndividualPeakScores(entry, spectrum, normalize, ions);
	}

	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Spectrum spectrum, boolean normalize, FragmentIon[] ions) {
		MassTolerance acquiredTolerance=parameters.getFragmentTolerance();
		MassTolerance libraryTolerance=parameters.getLibraryFragmentTolerance();
		
		double[] predictedMasses=entry.getMassArray();
		float[] predictedIntensities=entry.getIntensityArray();
		float[] correlation=entry.getCorrelationArray();
		
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();
		
		ArrayList<PeakScores> scoredPeaks=new ArrayList<PeakScores>();
		for (FragmentIon targetIon : ions) {
			double target=targetIon.mass;
			
			int[] predictedIndicies=libraryTolerance.getIndicies(predictedMasses, target);
			float predictedIntensity=0.0f;
			float maxCorrelation=0.01f;
			for (int i=0; i<predictedIndicies.length; i++) {
				if (predictedIntensity<predictedIntensities[predictedIndicies[i]]) {
					predictedIntensity=predictedIntensities[predictedIndicies[i]];
				}
				if (maxCorrelation<correlation[predictedIndicies[i]]) {
					maxCorrelation=correlation[predictedIndicies[i]];
				}
			}
			
			if (predictedIntensity>0) {
				int[] indicies=acquiredTolerance.getIndicies(acquiredMasses, target);
				float intensity=0.0f;
				float bestPeakIntensity=0.0f;
				float deltaMass=0.0f;
				for (int j=0; j<indicies.length; j++) {
					intensity+=acquiredIntensities[indicies[j]];
					
					if (acquiredIntensities[indicies[j]]>bestPeakIntensity) {
						bestPeakIntensity=acquiredIntensities[indicies[j]];

						deltaMass=(float)acquiredTolerance.getDeltaScore(target, acquiredMasses[indicies[j]]);
					}
				}
				float peakScore=predictedIntensity*intensity*maxCorrelation;
				if (intensity>0.0f) {
					scoredPeaks.add(new PeakScores(peakScore, targetIon, deltaMass));
				} else {
					scoredPeaks.add(null);
				}
			}
		}
		return scoredPeaks.toArray(new PeakScores[scoredPeaks.size()]);
	}
}
