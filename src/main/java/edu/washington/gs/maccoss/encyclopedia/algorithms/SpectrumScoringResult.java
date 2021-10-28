package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.AbstractRetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TFloatArrayList;

public class SpectrumScoringResult extends AbstractScoringResult {
	
	private final FragmentScan msms;
	private final ArrayList<Pair<ScoredObject<PeptidePrecursor>, float[]>> goodPeptides=new ArrayList<>(); // PeptidePrecursors are all LibraryEntries
	private XYTraceInterface trace=null;
	
	public SpectrumScoringResult(FragmentScan msms) {
		this.msms=msms;
	}
	
	public SpectrumScoringResult rescore(RetentionTimeAlignmentInterface filter) {
		SpectrumScoringResult newResult=new RescoredSpectrumScoringResult(msms);
		newResult.setTrace(trace);
		
		boolean anyFoundWithRTFilter=false;
		
		boolean bestSet=false;
		float bestScore=0;
		float[] bestScores=null;
		LibraryEntry bestPeptide=null;
		float bestActualRT=0f;
		float bestModelRT=0f;
		
		/* assumes sorted in order, first is best scoring! */
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
			
			if (!bestSet) {
				bestSet=true;
				bestScore=score;
				float deltaRT=Math.abs(filter.getDelta(actualRT, modelRT));
				float[] scoresWithRT=General.concatenate(scores, deltaRT);
				bestScores=scoresWithRT;
				bestPeptide=peptide;
				bestActualRT=actualRT;
				bestModelRT=modelRT;
			}
			
			boolean passes=filter.getProbabilityFitsModel(actualRT, modelRT)>=AbstractRetentionTimeFilter.rejectionPValue;
			float scorePercentage=score/bestScore;
			if (passes&&scorePercentage>0.9) {
				// next best score has to be close!
				float deltaRT=Math.abs(filter.getDelta(actualRT, modelRT));
				float[] scoresWithRT=General.concatenate(scores, deltaRT);
				newResult.addPeptide(score, scoresWithRT, peptide);
				anyFoundWithRTFilter=true;
			} 
		}
		
		// if nothing passes the RT filter then use the top match
		if (!anyFoundWithRTFilter) {
			// Spectrum-centric searching can just drop PSMs that don't fit the RT expectations, so could relax a little but require a match. So far this doesn't help
			//boolean passes=filter.getProbabilityFitsModel(bestActualRT, bestModelRT)>=AbstractRetentionTimeFilter.rejectionPValue;
			//if (passes) {
			newResult.addPeptide(bestScore, bestScores, bestPeptide);
			//}
		}
		
		return newResult;
	}
	
	public int size() {
		return goodPeptides.size();
	}
	
	public FragmentScan getMSMS() {
		return msms;
	}

	public void addPeptide(float score, float[] auxScoreArray, LibraryEntry peptide) {
		goodPeptides.add(new Pair<ScoredObject<PeptidePrecursor>, float[]>(new ScoredObject<PeptidePrecursor>(score, peptide), auxScoreArray));
	}
	
	@Override
	public void addStripe(float score, float[] auxScoreArray, FragmentScan stripe) {
		throw new EncyclopediaException("Unexpected addStripe in SpectrumScoringResult. You can only addPeptide to a SpectrumScoringResult (DDA)");	
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
		for (Pair<ScoredObject<PeptidePrecursor>, float[]> pair : getGoodPeptides()) {
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
		return goodPeptides.size()>0;
	}

	/**
	 * hack, can only report one good MSMS for a spectrum scoring result (since there's only one spectrum)
	 */
	public Pair<ScoredObject<FragmentScan>, float[]> getScoredMSMS() {
		Pair<ScoredObject<PeptidePrecursor>, float[]> pair=getScoredPeptide();
		return new Pair<ScoredObject<FragmentScan>, float[]>(new ScoredObject<FragmentScan>(pair.x.x, msms), pair.y);
	}
	
	/**
	 * hack, can only report one good MSMS for a spectrum scoring result (since there's only one spectrum)
	 */
	public ArrayList<Pair<ScoredObject<FragmentScan>, float[]>> getGoodMSMSCandidates() {
		ArrayList<Pair<ScoredObject<FragmentScan>, float[]>> list=new ArrayList<>();
		list.add(getScoredMSMS());
		return list;
	}
	
	/**
	 * hack, can only report one good MSMS for a spectrum scoring result (since there's only one spectrum)
	 */
	@Override
	public LibraryEntry getEntry() {
		PeptidePrecursor peptide = getScoredPeptide().x.y;
		if (peptide instanceof LibraryEntry) {
			return (LibraryEntry)peptide;
		} else {
			throw new EncyclopediaException("Expect LibraryEntry but got basic PeptidePrecursor in SpectrumScoringResult (DDA)");
		}
	}
	
	public ArrayList<Pair<ScoredObject<PeptidePrecursor>, float[]>> getGoodPeptides() {
		return goodPeptides;
	}
	
	public Pair<ScoredObject<PeptidePrecursor>, float[]> getScoredPeptide() {
		float bestScore=-Float.MAX_VALUE;
		Pair<ScoredObject<PeptidePrecursor>, float[]> bestPair=null;
		
		for (Pair<ScoredObject<PeptidePrecursor>, float[]> pair : goodPeptides) {
			if (pair.x.x>bestScore) {
				bestScore=pair.x.x;
				bestPair=pair;
			}
		}
		return bestPair;
	}
}
