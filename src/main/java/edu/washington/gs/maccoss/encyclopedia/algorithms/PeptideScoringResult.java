package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.AbstractRetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

public class PeptideScoringResult {
	public static final PeptideScoringResult POISON_RESULT=new PeptideScoringResult(null);
	
	private final LibraryEntry entry;
	private final ArrayList<Pair<ScoredObject<FragmentScan>, float[]>> goodStripes=new ArrayList<Pair<ScoredObject<FragmentScan>, float[]>>();
	private XYTraceInterface trace=null;
	
	public PeptideScoringResult(LibraryEntry entry) {
		this.entry=entry;
	}
	
	public PeptideScoringResult rescore(RetentionTimeAlignmentInterface filter) {
		PeptideScoringResult newResult=new RescoredPeptideScoringResult(entry);
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
		float bestScore=-Float.MAX_VALUE;
		for (Pair<ScoredObject<FragmentScan>, float[]> pair : goodStripes) {
			if (pair.x.x>bestScore) {
				bestScore=pair.x.x;
			}
		}
		return bestScore;
	}
	
	public void setTrace(XYTraceInterface trace) {
		this.trace=trace;
	}
	
	public XYTraceInterface getTrace() {
		return trace;
	}
	
	public ArrayList<Pair<ScoredObject<FragmentScan>, float[]>> getGoodStripes() {
		return goodStripes;
	}
}
