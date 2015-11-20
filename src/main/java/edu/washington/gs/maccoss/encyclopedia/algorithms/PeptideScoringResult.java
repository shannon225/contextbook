package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

public class PeptideScoringResult {
	private final ArrayList<ScoredObject<Stripe>> goodStripes=new ArrayList<ScoredObject<Stripe>>();
	private XYTrace trace=null;
	
	public PeptideScoringResult() {
	}

	public void addStripe(float score, Stripe stripe) {
		goodStripes.add(new ScoredObject<Stripe>(score, stripe));
	}
	
	public void setTrace(XYTrace trace) {
		this.trace=trace;
	}
	
	public XYTrace getTrace() {
		return trace;
	}
	
	public ArrayList<ScoredObject<Stripe>> getGoodStripes() {
		return goodStripes;
	}
}
