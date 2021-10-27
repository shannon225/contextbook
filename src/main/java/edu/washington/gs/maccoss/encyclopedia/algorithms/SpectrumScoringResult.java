package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.AbstractRetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

public class SpectrumScoringResult {
	public static final SpectrumScoringResult POISON_RESULT=new SpectrumScoringResult(null);
	
	private final FragmentScan msms;
	private final ArrayList<Pair<ScoredObject<PeptidePrecursor>, float[]>> goodPeptides=new ArrayList<>();
	private XYTraceInterface trace=null;
	
	public SpectrumScoringResult(FragmentScan msms) {
		this.msms=msms;
	}
	
	public SpectrumScoringResult rescore(RetentionTimeAlignmentInterface filter) {
		SpectrumScoringResult newResult=new RescoredSpectrumScoringResult(msms);
		newResult.setTrace(trace);
		
		boolean anyFoundWithRTFilter=false;
		boolean bestSet=false;
		float bestScore=0.0f;
		float[] bestScores=null;
		PeptidePrecursor bestPeptide=null;
		
		for (Pair<ScoredObject<PeptidePrecursor>, float[]> pair : goodPeptides) {
			float score=pair.x.x;
			PeptidePrecursor precursor=pair.x.y;
			if (!(precursor instanceof LibraryEntry)) {
				continue;
			}
			LibraryEntry peptide=(LibraryEntry)precursor;
			float[] scores=pair.y;
			float modelRT=peptide.getScanStartTime()/60f;
			float actualRT=msms.getScanStartTime()/60f;
			boolean passes=filter.getProbabilityFitsModel(actualRT, modelRT)>=AbstractRetentionTimeFilter.rejectionPValue;
			if (passes) {
				float deltaRT=Math.abs(actualRT-filter.getYValue(modelRT));
				float[] scoresWithRT=General.concatenate(scores, deltaRT);
				newResult.addPeptide(score, scoresWithRT, peptide);
				anyFoundWithRTFilter=true;
			} else if (!bestSet) {
				bestSet=true;
				bestScore=score;
				float deltaRT=Math.abs(actualRT-filter.getYValue(modelRT));
				float[] scoresWithRT=General.concatenate(scores, deltaRT);
				bestScores=scoresWithRT;
				bestPeptide=peptide;
			}
		}
		
		// if nothing passes the RT filter then use the top match
		if (!anyFoundWithRTFilter) {
			newResult.addPeptide(bestScore, bestScores, bestPeptide);
		}
		
		return newResult;
	}
	
	public int size() {
		return goodPeptides.size();
	}
	
	public FragmentScan getMSMS() {
		return msms;
	}

	public void addPeptide(float score, float[] auxScoreArray, PeptidePrecursor Peptide) {
		goodPeptides.add(new Pair<ScoredObject<PeptidePrecursor>, float[]>(new ScoredObject<PeptidePrecursor>(score, Peptide), auxScoreArray));
	}
	
	public float getBestScore() {
		float bestScore=-Float.MAX_VALUE;
		for (Pair<ScoredObject<PeptidePrecursor>, float[]> pair : goodPeptides) {
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
	
	public ArrayList<Pair<ScoredObject<PeptidePrecursor>, float[]>> getGoodPeptides() {
		return goodPeptides;
	}
}
