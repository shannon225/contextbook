package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneAuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;

public class XCorDIAOneScorer implements PSMScorer {
	private final EncyclopediaOneAuxillaryPSMScorer auxScorer;

	public XCorDIAOneScorer(SearchParameters parameters, LibraryBackgroundInterface background) {
		auxScorer=new EncyclopediaOneAuxillaryPSMScorer(parameters, background, false);
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
	public float score(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return score(entry, spectrum);
	}

	public float score(LibraryEntry entry, Spectrum spectrum) {
		if (entry instanceof XCorrLibraryEntry&&spectrum instanceof XCorrStripe) {
			return ((XCorrLibraryEntry) entry).score(((XCorrStripe)spectrum).getXcorrSpectrum());
		}
		throw new EncyclopediaException("XCorDIAOneScorer can only score XCorr-based data types");
	}
}
