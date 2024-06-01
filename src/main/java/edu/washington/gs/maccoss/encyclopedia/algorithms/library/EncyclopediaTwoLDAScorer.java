package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.TargeteDecoyPSMFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;

public class EncyclopediaTwoLDAScorer implements EncyclopediaScorer {
	private final EncyclopediaTwoScorer mainScorer;
	private final Optional<Pair<MProphetResult, TargeteDecoyPSMFilter>> ldaScorer;
	private final int[] scoreIndexConvertingArray;
	
	public EncyclopediaTwoLDAScorer(EncyclopediaTwoScorer mainScorer, Optional<Pair<MProphetResult, TargeteDecoyPSMFilter>> ldaScorer) {
		this.mainScorer = mainScorer;
		this.ldaScorer = ldaScorer;
		
		if (!ldaScorer.isPresent()) {
			scoreIndexConvertingArray=null;
		} else {
			ArrayList<String> featureNames=ldaScorer.get().getX().getFeatureNames();
			scoreIndexConvertingArray=new int[featureNames.size()];
			Arrays.fill(scoreIndexConvertingArray, -1);

			String[] auxScoreNames=EncyclopediaTwoAuxillaryPSMScorer.getScoreNames();
			for (int i = 0; i < scoreIndexConvertingArray.length; i++) {
				String scoreName=featureNames.get(i);
				for (int j = 0; j < auxScoreNames.length; j++) {
					if (scoreName.equals(auxScoreNames[j])) {
						scoreIndexConvertingArray[i]=j;
						break;
					}
				}
			}
		}
	}
	
	public Optional<Pair<MProphetResult, TargeteDecoyPSMFilter>> getLDAScorer() {
		return ldaScorer;
	}
	public int[] getScoreIndexConvertingArray() {
		return scoreIndexConvertingArray;
	}

	public static String getPrimaryScoreName() {
		return "primary";
	}

	@Override
	public EncyclopediaTwoAuxillaryPSMScorer getAuxScorer() {
		return mainScorer.getAuxScorer();
	}

	@Override
	public float score(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return mainScorer.score(entry, spectrum, predictedIsotopeDistribution, precursors);
	}
	
	@Override
	public float[] auxScore(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return mainScorer.getAuxScorer().score(entry, spectrum, predictedIsotopeDistribution, precursors);
	}
	@Override
	public String[] getAuxScoreNames(LibraryEntry entry) {
		return mainScorer.getAuxScorer().getScoreNames(entry);
	}

	@Override
	public float score(LibraryEntry entry, Spectrum spectrum) {
		PeakScores[] individualPeakScores=getIndividualPeakScores(entry, spectrum, true);
		return mainScorer.scoreIons(individualPeakScores);
	}

	@Override
	public int getParentDeltaMassIndex() {
		return mainScorer.getAuxScorer().getParentDeltaMassIndex();
	}

	@Override
	public int getFragmentDeltaMassIndex() {
		return mainScorer.getAuxScorer().getFragmentDeltaMassIndex();
	}

	@Override
	public float score(LibraryEntry entry, Spectrum spectrum, FragmentIon[] ions) {
		PeakScores[] individualPeakScores=getIndividualPeakScores(entry, spectrum, true, ions);
		return mainScorer.scoreIons(individualPeakScores);
	}

	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Spectrum spectrum, boolean normalize) {
		return mainScorer.getIndividualPeakScores(entry, spectrum, normalize);
	}

	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Spectrum spectrum, boolean normalize, FragmentIon[] ions) {
		return mainScorer.getIndividualPeakScores(entry, spectrum, normalize, ions);
	}
}
