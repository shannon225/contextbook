package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import gnu.trove.map.hash.TDoubleObjectHashMap;

public class PecanOneScoringFactory implements PecanScoringFactory {
	public static final String version="0.3.3";
	private final PecanSearchParameters parameters;
	private final File featureFile;

	public PecanOneScoringFactory(PecanSearchParameters parameters, File featureFile) {
		this.parameters=parameters;
		this.featureFile=featureFile;
	}
	
	public String getVersion() {
		return version;
	}
	
	@Override
	public PecanSearchParameters getParameters() {
		return parameters;
	}
	
	@Override
	public AbstractPecanFragmentationModel getFragmentationModel(FastaPeptideEntry peptide, AminoAcidConstants aaConstants) {
		return new PecanOneFragmentationModel(peptide, aaConstants);
	}

	@Override
	public PSMScorer getBackgroundScorer() {
		return new PecanRawScorer(parameters.getFragmentTolerance(), null);
	}

	@Override
	public PSMPeakScorer getPecanScorer() {
		return new PecanRawScorer(parameters.getFragmentTolerance(), new PecanAuxillaryScorer(parameters));
	}

	@Override
	public AbstractPecanScoringTask getScoringTask(PSMPeakScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, TDoubleObjectHashMap<XYPoint>[] background, PrecursorScanMap precursors,
			int scanAveragingMargin, BlockingQueue<PeptideScoringResult> resultsQueue) {
		return new PecanOneScoringTask(scorer, entries, stripes, background, precursors, scanAveragingMargin, resultsQueue, parameters);
	}

	@Override
	public PeptideScoringResultsConsumer getResultsConsumer(BlockingQueue<PeptideScoringResult> resultsQueue, StripeFileInterface diaFile) {
		return new PecanScoringResultsToTSVConsumer(featureFile, diaFile, resultsQueue, parameters.getNumberOfReportedPeaks());
	}
}
