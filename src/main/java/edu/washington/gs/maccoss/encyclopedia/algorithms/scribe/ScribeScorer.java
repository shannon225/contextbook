package edu.washington.gs.maccoss.encyclopedia.algorithms.scribe;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorrLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorrStripe;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;

public class ScribeScorer implements PSMScorer {
	private final ScribeAuxillaryPSMScorer auxScorer;

	public ScribeScorer(SearchParameters parameters) {
		auxScorer=new ScribeAuxillaryPSMScorer(parameters);
	}
	
	@Override
	public float[] auxScore(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return auxScorer.score(entry, spectrum, predictedIsotopeDistribution, precursors);
	}
	
	@Override
	public float score(LibraryEntry entry, Spectrum spectrum, FragmentIon[] ions) {
		return score(entry, spectrum);
	}

	@Override
	public String[] getAuxScoreNames(LibraryEntry entry) {
		return auxScorer.getScoreNames(entry);
	}

	@Override
	public int getParentDeltaMassIndex() {
		return auxScorer.getParentDeltaMassIndex();
	}

	@Override
	public int getFragmentDeltaMassIndex() {
		return auxScorer.getFragmentDeltaMassIndex();
	}

	@Override
	public float score(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return score(entry, spectrum);
	}

	public float score(LibraryEntry entry, Spectrum spectrum) {
		if (entry instanceof XCorrLibraryEntry&&spectrum instanceof XCorrStripe) {
			return ((XCorrLibraryEntry) entry).score(((XCorrStripe)spectrum).getXcorrSpectrum());
		}
		throw new EncyclopediaException("ScribeScorer can only score XCorr-based data types");
	}
}
