package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SpectrumUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.set.hash.TIntHashSet;

public class XCorDIAOneScoringTask extends AbstractLibraryScoringTask {
	private final float dutyCycle;
	
	public XCorDIAOneScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue,
			SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		this.dutyCycle=dutyCycle;
		
		assert(scorer instanceof XCorDIAOneScorer);
	}
	
	private static final int peaksKept=0;

	@Override
	protected Nothing process() {
		int movingAverageLength=Math.round(parameters.getExpectedPeakWidth()/dutyCycle);
		for (LibraryEntry entry : super.entries) {
			XCorrLibraryEntry xcordiaEntry;
			if (entry instanceof XCorrLibraryEntry) {
				xcordiaEntry=(XCorrLibraryEntry)entry;
			} else {
				xcordiaEntry=XCorrLibraryEntry.generateEntry(false, entry.getSource(), entry.getAccessions(), entry.getPrecursorCharge(), entry.getPeptideModSeq(), parameters);
				
			}
			xcordiaEntry.init();
			
			PeptideScoringResult result=new PeptideScoringResult(entry);
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq(), parameters.getAAConstants());
			
			float[] primary=new float[super.stripes.size()];
			float[] rts=new float[super.stripes.size()];
			for (int i=0; i<super.stripes.size(); i++) {
				Stripe stripe=super.stripes.get(i);
				XCorrStripe xcordiaStripe;
				if (stripe instanceof XCorrStripe) {
					xcordiaStripe=(XCorrStripe)stripe;
				} else {
					xcordiaStripe=new XCorrStripe(stripe, parameters);
				}
				primary[i]=scorer.score(xcordiaEntry, xcordiaStripe, predictedIsotopeDistribution, precursors);
				
				rts[i]=xcordiaStripe.getScanStartTime();
			}
			
			float[] averagePrimary=gaussianCenteredAverage(primary, movingAverageLength);

			TFloatFloatHashMap map=new TFloatFloatHashMap();
			ArrayList<ScoredIndex> goodStripes=new ArrayList<ScoredIndex>();
			for (int i=0; i<averagePrimary.length; i++) {
				goodStripes.add(new ScoredIndex(primary[i], i));
				map.put(i, primary[i]);
			}
			Collections.sort(goodStripes);

			EValueCalculator calculator=new EValueCalculator(map);

			TIntHashSet takenScans=new TIntHashSet();
			int identifiedPeaks=0;
			for (int i=goodStripes.size()-1; i>=0; i--) {
				float score=goodStripes.get(i).x;
				int index=goodStripes.get(i).y;
				if (takenScans.contains(index)) {
					continue;
					
				} else {
					Stripe stripe=super.stripes.get(index);
					float[] auxScoreArray=scorer.auxScore(entry, stripe, predictedIsotopeDistribution, precursors);
					float evalue=calculator.getNegLog10EValue(score);
					if (Float.isNaN(evalue)) {
						evalue=-1.0f;
					}
					result.addStripe(score, General.concatenate(auxScoreArray, score, evalue), stripe);
					
					// block out a 40 scan window
					int lowerWindow=index-2*movingAverageLength;
					int upperWindow=index+2*movingAverageLength;
					for (int j=lowerWindow; j<=upperWindow; j++) {
						takenScans.add(j);
					}
					
					if (identifiedPeaks>peaksKept) {
						// keep N+1 peaks
						break;
					}
					identifiedPeaks++;
				}
			}

			XYTrace trace=new XYTrace(rts, primary, GraphType.line, "XCorr");
			result.setTrace(trace);
			resultsQueue.add(result);
		}
		return Nothing.NOTHING;
	}

	
}
