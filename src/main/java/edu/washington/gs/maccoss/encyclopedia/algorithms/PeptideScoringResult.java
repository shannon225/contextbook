package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

public class PeptideScoringResult {
	private final ArrayList<Pair<ScoredObject<Stripe>, float[]>> goodStripes=new ArrayList<Pair<ScoredObject<Stripe>, float[]>>();
	private XYTrace trace=null;
	
	public PeptideScoringResult() {
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
