package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

public class PeptideScoringResult {
	public static final PeptideScoringResult POISON_RESULT=new PeptideScoringResult(null);
	
	private final LibraryEntry entry;
	private final ArrayList<Pair<ScoredObject<Stripe>, float[]>> goodStripes=new ArrayList<Pair<ScoredObject<Stripe>, float[]>>();
	private XYTrace trace=null;
	
	public PeptideScoringResult(LibraryEntry entry) {
		this.entry=entry;
	}
	
	public PeptideScoringResult rescore(RetentionTimeFilter filter) {
		PeptideScoringResult newResult=new RescoredPeptideScoringResult(entry);
		newResult.setTrace(trace);
		
		for (Pair<ScoredObject<Stripe>, float[]> pair : goodStripes) {
			float score=pair.x.x;
			Stripe stripe=pair.x.y;
			float[] scores=pair.y;
			float actualRT=stripe.getScanStartTime()/60f;
			float modelRT=entry.getRetentionTime();
			boolean passes=filter.getProbabilityFitsModel(actualRT, modelRT)>=0.05f;
			if (passes) {
				float deltaRT=Math.abs(actualRT-filter.getYValue(modelRT));
				float[] scoresWithRT=General.concatenate(scores, deltaRT);
				newResult.addStripe(score, scoresWithRT, stripe);
			}
		}
		
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
