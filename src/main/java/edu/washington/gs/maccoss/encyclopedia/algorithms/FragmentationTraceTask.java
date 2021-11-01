package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanRawScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;

public class FragmentationTraceTask extends PeptideScoringTask {
	public static final byte PLOT_INTENSITIES=0;
	public static final byte PLOT_SCORES=1;
	public static final byte PLOT_DELTA_MASSES=2;
	
	private final byte plottingMethod;

	public FragmentationTraceTask(PecanRawScorer scorer, byte plottingMethod, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, PrecursorScanMap precursors, AminoAcidConstants aaConstants) {
		super(scorer, entries, stripes, precursors, aaConstants);
		this.plottingMethod=plottingMethod;
	}
	private PecanRawScorer getScorer() {
		return (PecanRawScorer)super.scorer;
	}

	@Override
	protected HashMap<LibraryEntry, AbstractScoringResult> process() {
		HashMap<LibraryEntry, AbstractScoringResult> map=new HashMap<LibraryEntry, AbstractScoringResult>();

		for (LibraryEntry entry : super.entries) {
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq(), aaConstants);
			String[] scoreNames=getScorer().getAuxScoreNames(entry);
			
			@SuppressWarnings("unchecked")
			ArrayList<XYPoint>[] dataPoints=new ArrayList[scoreNames.length];
			for (int i=0; i<dataPoints.length; i++) {
				dataPoints[i]=new ArrayList<XYPoint>();
			}
			for (FragmentScan stripe : super.stripes) {
				float rt=stripe.getScanStartTime()/60.0f;
				
				if (plottingMethod==PLOT_INTENSITIES) {
					float[] fragmentScores=getScorer().auxScore(entry, stripe, predictedIsotopeDistribution, precursors);
					for (int i=0; i<fragmentScores.length; i++) {
						dataPoints[i].add(new XYPoint(rt, fragmentScores[i]));
					}
				} else {
					PeakScores[] scores=getScorer().getIndividualPeakScores(entry, stripe, false);
					
					if (plottingMethod==PLOT_DELTA_MASSES) {
						for (int i=0; i<scores.length; i++) {
							if (scores[i]!=null) {
								dataPoints[i].add(new XYPoint(rt, scores[i].getDeltaMass()));
							} else {
								dataPoints[i].add(new XYPoint(rt, Double.NaN));
							}
						}
					} else if (plottingMethod==PLOT_SCORES) {
						for (int i=0; i<scores.length; i++) {
							if (scores[i]!=null) {
								dataPoints[i].add(new XYPoint(rt, scores[i].getScore()));
							} else {
								dataPoints[i].add(new XYPoint(rt, Double.NaN));
							}
						}
					}
				}
			}
			
			FragmentationScoringResult result=new FragmentationScoringResult(entry);
			for (int i=0; i<dataPoints.length; i++) {
				XYTrace trace=new XYTrace(dataPoints[i], GraphType.line, scoreNames[i]);
				result.addFragmentationTrace(trace);
			}

			map.put(entry, result);
		}
		return map;
	}
}
