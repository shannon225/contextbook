package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RunningMedianWarper;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

public class PeptideScoringResult {
	public static final PeptideScoringResult POISON_RESULT=new PeptideScoringResult(null);
	
	private final LibraryEntry entry;
	private final ArrayList<Pair<ScoredObject<Stripe>, float[]>> goodStripes=new ArrayList<Pair<ScoredObject<Stripe>, float[]>>();
	private XYTrace trace=null;
	
	public PeptideScoringResult(LibraryEntry entry) {
		this.entry=entry;
	}
	
	public PeptideScoringResult rescore(Pair<RunningMedianWarper, LinearDiscriminantAnalysis> rescoringModel) {
		PeptideScoringResult newResult=new RescoredPeptideScoringResult(entry);
		newResult.setTrace(trace);
		
		/*float bestScore=-Float.MAX_VALUE;
		float[] bestAuxs=null;
		Stripe bestStripe=null;
		*/
		
		for (Pair<ScoredObject<Stripe>, float[]> pair : goodStripes) {
			Stripe stripe=pair.x.y;
			float[] scores=pair.y;
			float entryTime=rescoringModel.x.getYValue(entry.getRetentionTime());
			float deltaRT=stripe.getScanStartTime()/60f-entryTime;
			float[] scoresWithRT=General.concatenate(scores, deltaRT);
			float newScore=rescoringModel.y.getScore(scoresWithRT);

			scoresWithRT=General.concatenate(scores, deltaRT, newScore);
			newResult.addStripe(newScore, scoresWithRT, stripe);
			/*if (newScore>bestScore) {
				bestScore=newScore;
				bestAuxs=scoresWithRT;
				bestStripe=stripe;
			}*/
		}
		/*if (bestStripe!=null) {
			newResult.addStripe(bestScore, bestAuxs, bestStripe);
		}*/
		
		return newResult;
	}
	
	public int size() {
		return goodStripes.size();
	}
	
	public LibraryEntry getEntry() {
		return entry;
	}

	public void addStripe(float score, float[] auxScoreArray, Stripe stripe) {
		goodStripes.add(new Pair<ScoredObject<Stripe>, float[]>(new ScoredObject<Stripe>(score, stripe), auxScoreArray));
	}
	
	public void setTrace(XYTrace trace) {
		this.trace=trace;
	}
	
	public XYTrace getTrace() {
		return trace;
	}
	
	public ArrayList<Pair<ScoredObject<Stripe>, float[]>> getGoodStripes() {
		return goodStripes;
	}
}
