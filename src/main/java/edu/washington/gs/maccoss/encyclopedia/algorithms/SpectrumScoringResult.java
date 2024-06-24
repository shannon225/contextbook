package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.ScoredPSMFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.ScoredPSMFilterInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;

public class SpectrumScoringResult extends AbstractScoringResult {
	
	private final FragmentScan msms;
	private final ArrayList<ScoredPSM> goodPeptides=new ArrayList<>(); // PeptidePrecursors are all LibraryEntries
	private XYTraceInterface trace=null;
	
	public SpectrumScoringResult(FragmentScan msms) {
		this.msms=msms;
	}
	
	public SpectrumScoringResult rescore(ScoredPSMFilterInterface filter) {
		SpectrumScoringResult newResult;
		if (filter instanceof ScoredPSMFilter) {
			newResult=new RecalibratedSpectrumScoringResult(msms);
		} else {
			newResult=new RescoredSpectrumScoringResult(msms);
		}
		newResult.setTrace(trace);
		
		boolean anyFoundWithRTFilter=false;
		
		if (goodPeptides.size()==0) return newResult;
		
		/* assumes sorted in order, first is best scoring! */
		ScoredPSM startingBest=goodPeptides.get(0);
				
		for (ScoredPSM pair : goodPeptides) {
			// next best score has to be 90% close!
			if (filter.passesFilter(pair)&&pair.getPrimaryScore()/startingBest.getPrimaryScore()>0.9f) {
				float[] scoresWithRT=General.concatenate(pair.getAuxScores(), filter.getAdditionalScores(pair));
				newResult.addPeptide(pair.getPrimaryScore(), scoresWithRT, pair.getDeltaPrecursorMass(), pair.getDeltaFragmentMass(), pair.getLibraryEntry());
				anyFoundWithRTFilter=true;
			}
		}
		
		// if nothing passes the RT filter then use the top match
		if (!anyFoundWithRTFilter) {
			float[] scoresWithRT=General.concatenate(startingBest.getAuxScores(), filter.getAdditionalScores(startingBest));
			newResult.addPeptide(startingBest.getPrimaryScore(), scoresWithRT, startingBest.getDeltaPrecursorMass(), startingBest.getDeltaFragmentMass(), startingBest.getLibraryEntry());
		}
		
		return newResult;
	}
	
	public int size() {
		return goodPeptides.size();
	}
	
	public FragmentScan getMSMS() {
		return msms;
	}

	public void addPeptide(float score, float[] auxScoreArray, float deltaPrecursorMass, float deltaFragmentMass, LibraryEntry peptide) {
		goodPeptides.add(new ScoredPSM(peptide, msms, score, auxScoreArray, deltaPrecursorMass, deltaFragmentMass));
	}
	
	@Override
	public void addStripe(float score, float[] auxScoreArray, float deltaPrecursorMass, float deltaFragmentMass, FragmentScan stripe) {
		throw new EncyclopediaException("Unexpected addStripe in SpectrumScoringResult. You can only addPeptide to a SpectrumScoringResult (DDA)");	
	}
	
	public void trim(int trimToN) {
		Collections.sort(goodPeptides);
		Collections.reverse(goodPeptides);
		while (goodPeptides.size()>trimToN) {
			goodPeptides.remove(goodPeptides.size()-1);
		}
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
		for (ScoredPSM pair : getGoodPeptides()) {
			scores.add(pair.getPrimaryScore());
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
	public ScoredPSM getScoredMSMS() {
		return getScoredPeptide();
	}
	
	/**
	 * hack, can only report one good MSMS for a spectrum scoring result (since there's only one spectrum)
	 */
	public ArrayList<ScoredPSM> getGoodMSMSCandidates() {
		ArrayList<ScoredPSM> list=new ArrayList<>();
		list.add(getScoredMSMS());
		return list;
	}
	
	/**
	 * hack, can only report one good MSMS for a spectrum scoring result (since there's only one spectrum)
	 */
	@Override
	public LibraryEntry getEntry() {
		return getScoredPeptide().getLibraryEntry();
	}
	
	public ArrayList<ScoredPSM> getGoodPeptides() {
		return goodPeptides;
	}
	
	public ScoredPSM getScoredPeptide() {
		float bestScore=-Float.MAX_VALUE;
		ScoredPSM bestPair=null;
		
		for (ScoredPSM pair : goodPeptides) {
			if (pair.getPrimaryScore()>bestScore) {
				bestScore=pair.getPrimaryScore();
				bestPair=pair;
			}
		}
		return bestPair;
	}
}
