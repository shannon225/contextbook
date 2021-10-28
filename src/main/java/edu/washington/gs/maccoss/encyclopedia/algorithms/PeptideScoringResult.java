package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.AbstractRetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TFloatArrayList;

public class PeptideScoringResult extends AbstractScoringResult {
	
	private final LibraryEntry entry;
	private final ArrayList<Pair<ScoredObject<FragmentScan>, float[]>> goodStripes=new ArrayList<Pair<ScoredObject<FragmentScan>, float[]>>();
	private XYTraceInterface trace=null;
	
	public PeptideScoringResult(LibraryEntry entry) {
		this.entry=entry;
	}
	
	public AbstractScoringResult rescore(RetentionTimeAlignmentInterface filter) {
		AbstractScoringResult newResult=new RescoredPeptideScoringResult(entry);
		newResult.setTrace(trace);
		
		boolean anyFoundWithRTFilter=false;
		boolean bestSet=false;
		float bestScore=0.0f;
		float[] bestScores=null;
		FragmentScan bestStripe=null;
		
		for (Pair<ScoredObject<FragmentScan>, float[]> pair : goodStripes) {
			float score=pair.x.x;
			FragmentScan stripe=pair.x.y;
			float[] scores=pair.y;
			float actualRT=stripe.getScanStartTime()/60f;
			float modelRT=entry.getRetentionTime()/60f;
			boolean passes=filter.getProbabilityFitsModel(actualRT, modelRT)>=AbstractRetentionTimeFilter.rejectionPValue;
			if (passes) {
				float deltaRT=Math.abs(actualRT-filter.getYValue(modelRT));
				float[] scoresWithRT=General.concatenate(scores, deltaRT);
				newResult.addStripe(score, scoresWithRT, stripe);
				anyFoundWithRTFilter=true;
			} else if (!bestSet) {
				bestSet=true;
				bestScore=score;
				float deltaRT=Math.abs(actualRT-filter.getYValue(modelRT));
				float[] scoresWithRT=General.concatenate(scores, deltaRT);
				bestScores=scoresWithRT;
				bestStripe=stripe;
			}
		}
		
		// if nothing passes the RT filter then use the top match
		if (!anyFoundWithRTFilter) {
			newResult.addStripe(bestScore, bestScores, bestStripe);
		}
		
		return newResult;
	}
	
	public int size() {
		return goodStripes.size();
	}
	
	public LibraryEntry getEntry() {
		return entry;
	}

	public void addStripe(float score, float[] auxScoreArray, FragmentScan stripe) {
		goodStripes.add(new Pair<ScoredObject<FragmentScan>, float[]>(new ScoredObject<FragmentScan>(score, stripe), auxScoreArray));
	}
	
	public float getBestScore() {
		float[] scores=getSortedScores();
		if (scores.length>0) {
			return scores[scores.length-1];
		}
		return 0.0f;
	}
	public float getSecondBestScore() {
		float[] scores=getSortedScores();
		if (scores.length>1) {
			return scores[scores.length-2];
		}
		return 0.0f;
	}
	
	private float[] getSortedScores() {
		TFloatArrayList scores=new TFloatArrayList();
		for (Pair<ScoredObject<FragmentScan>, float[]> pair : goodStripes) {
			scores.add(pair.x.x);
		}
		float[] sorted=scores.toArray();
		Arrays.sort(sorted);
		return sorted;
	}
	
	public void setTrace(XYTraceInterface trace) {
		this.trace=trace;
	}
	
	public XYTraceInterface getTrace() {
		return trace;
	}
	
	public boolean hasScoredResults() {
		return goodStripes.size()>0;
	}
	
	public Pair<ScoredObject<FragmentScan>, float[]> getScoredMSMS() {
		float bestScore=-Float.MAX_VALUE;
		Pair<ScoredObject<FragmentScan>, float[]> bestPair=null;
		
		for (Pair<ScoredObject<FragmentScan>, float[]> pair : goodStripes) {
			if (pair.x.x>bestScore) {
				bestScore=pair.x.x;
				bestPair=pair;
			}
		}
		return bestPair;
	}
	
	public ArrayList<Pair<ScoredObject<FragmentScan>, float[]>> getGoodMSMSCandidates() {
		return goodStripes;
	}
}

