package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class EncyclopediaOneScorer implements EncyclopediaScorer {
	private final SearchParameters parameters;
	private final EncyclopediaOneAuxillaryPSMScorer auxScorer;

	public EncyclopediaOneScorer(SearchParameters parameters, LibraryBackgroundInterface background) {
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

	public PeakScores[] altgetIndividualPeakScores(LibraryEntry entry, Spectrum spectrum, boolean normalize, FragmentIon[] ions) {
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

	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Spectrum spectrum, boolean normalize, FragmentIon[] ions) {
		MassTolerance acquiredTolerance=parameters.getFragmentTolerance();
		MassTolerance libraryTolerance=parameters.getLibraryFragmentTolerance();
		
		double[] predictedMasses=entry.getMassArray();
		float[] predictedIntensities=entry.getIntensityArray();
		float[] correlation=entry.getCorrelationArray();
		
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();

		final List<FragmentIon> uniqueFragments;
		{
			final List<FragmentIon> work = Lists.newArrayList();

			Arrays.stream(ions)
					.sorted(Ordering.natural().onResultOf(FragmentIon::getType)) // sort by declaration order of fragment types -- this is roughly by priority
					.forEach(ion -> {
						boolean exists = false;
						for (FragmentIon existing : work) {
							if (acquiredTolerance.equals(ion.mass, existing.mass)) {
								Logger.logLine("Skipping fragment " + ion.toString() + " in favor of " + existing.toString() + " from " + entry.getPeptideModSeq() + "+" + entry.getPrecursorCharge());
								exists = true;
								break;
							}
						}

						if (!exists) {
							work.add(ion);
						}
					});

			uniqueFragments = ImmutableList.copyOf(work);
		}
		ions = uniqueFragments.toArray(new FragmentIon[uniqueFragments.size()]); // slightly naughty to overwrite parameter, but this is fine
		Arrays.sort(ions); // sort by natural ordering of ions

		ArrayList<PeakScores> scoredPeaks=new ArrayList<PeakScores>();
		
		int predictedIndex=0;
		int acquiredIndex=0;
		for (int i=0; i<ions.length; i++) {
			float predictedIntensity=0.0f;
			float maxCorrelation=0.01f;
			
			for (int j=predictedIndex; j<predictedMasses.length; j++) {
				int compare=libraryTolerance.compareTo(ions[i].mass, predictedMasses[j]);
				if (compare<0) {
					predictedIndex=j;
					break;
				} else if (compare==0) {
					if (predictedIntensity<predictedIntensities[j]) {
						predictedIntensity=predictedIntensities[j];
					}
					if (maxCorrelation<correlation[j]) {
						maxCorrelation=correlation[j];
					}
				}
			}

			float acquiredIntensity=0.0f;
			float deltaMass=0.0f;
			if (predictedIntensity>0.0f) {
				for (int j=acquiredIndex; j<acquiredMasses.length; j++) {
					int compare=libraryTolerance.compareTo(ions[i].mass, acquiredMasses[j]);
					if (compare<0) {
						acquiredIndex=j;
						break;
					} else if (compare==0) {
						if (acquiredIntensity<acquiredIntensities[j]) {
							acquiredIntensity=acquiredIntensities[j];
							deltaMass=(float)acquiredTolerance.getDeltaScore(ions[i].mass, acquiredMasses[j]);
						}
					}
				}
			}

			if (predictedIntensity>0.0f&&acquiredIntensity>0.0f) {
				float peakScore=predictedIntensity*acquiredIntensity*maxCorrelation;
				scoredPeaks.add(new PeakScores(peakScore, ions[i], deltaMass));
			} else {
				scoredPeaks.add(null);
			}
		}
		return scoredPeaks.toArray(new PeakScores[scoredPeaks.size()]);
	}
}
