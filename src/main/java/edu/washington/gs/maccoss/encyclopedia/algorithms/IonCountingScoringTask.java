package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.math.BackgroundSubtractionFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class IonCountingScoringTask extends AbstractLibraryScoringTask {
	private final float dutyCycle;
	private final int peaksKept=1;
	
	public IonCountingScoringTask(PSMPeakScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<AbstractScoringResult> resultsQueue,
			SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		this.dutyCycle=dutyCycle;
	}

	@Override
	protected Nothing process() {
		MassTolerance tolerance=parameters.getFragmentTolerance();
		int movingAverageLength=Math.round(parameters.getExpectedPeakWidth()/dutyCycle);
		for (LibraryEntry entry : super.entries) {
			AbstractScoringResult result=new PeptideScoringResult(entry);

			TFloatFloatHashMap scoreMap=new TFloatFloatHashMap();
			float[] primary=new float[super.stripes.size()];
			for (int i=0; i<super.stripes.size(); i++) {
				FragmentScan stripe=super.stripes.get(i);
				double[] masses=stripe.getMassArray();
				
				for (double mass : entry.getMassArray()) {
					int[] indicies=tolerance.getIndices(masses, mass);
					if (indicies.length>0) primary[i]+=1.0f;
				}

				scoreMap.put(stripe.getScanStartTime(), primary[i]);
			}
			
			float[] averagePrimary=BackgroundSubtractionFilter.gaussianCenteredAverage(primary, movingAverageLength);

			TFloatFloatHashMap map=new TFloatFloatHashMap();
			ArrayList<ScoredIndex> goodStripes=new ArrayList<ScoredIndex>();
			for (int i=0; i<averagePrimary.length; i++) {
				goodStripes.add(new ScoredIndex(primary[i], i));
				map.put(i, primary[i]);
			}
			Collections.sort(goodStripes);

			int identifiedPeaks=0;
			for (int i=goodStripes.size()-1; i>=0; i--) {
				float score=goodStripes.get(i).x;
				int index=goodStripes.get(i).y;

				FragmentScan stripe=super.stripes.get(index);
				result.addStripe(score, new float[0], 0.0f, 0.0f, stripe);
				
				
				if (identifiedPeaks>peaksKept) {
					// keep N+1 peaks
					break;
				}
				identifiedPeaks++;
			}
			
			result.setTrace(new XYTrace(scoreMap, GraphType.line, entry.getPeptideModSeq()));
			
			resultsQueue.add(result);
		}
		return Nothing.NOTHING;
	}

	
}
