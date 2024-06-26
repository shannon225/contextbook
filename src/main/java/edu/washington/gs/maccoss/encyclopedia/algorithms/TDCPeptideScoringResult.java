package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.ScoredPSMFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.ScoredPSMFilterInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.TargeteDecoyPSMFilter;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;

public class TDCPeptideScoringResult extends AbstractScoringResult {
	
	private final ArrayList<ScoredPSM> goodStripes=new ArrayList<ScoredPSM>();
	
	public TDCPeptideScoringResult() {
	}
	
	public AbstractScoringResult rescore(ScoredPSMFilterInterface filter) {
		
		if (goodStripes.size()==0) return null;
		
		LibraryEntry entry=getEntry();
		AbstractScoringResult newResult;
		if (filter instanceof ScoredPSMFilter||filter instanceof TargeteDecoyPSMFilter) {
			newResult=new RecalibratedPeptideScoringResult(entry);
		} else {
			newResult=new RescoredPeptideScoringResult(entry);
		}
		
		boolean anyFoundWithRTFilter=false;
		
		/* assumes sorted in order, first is best scoring! */
		ScoredPSM startingBest=goodStripes.get(0);
				
		for (ScoredPSM pair : goodStripes) {
			// next best score has to be 90% close!
			if (filter.passesFilter(pair)&&pair.getPrimaryScore()/startingBest.getPrimaryScore()>0.9f) {
				float[] scoresWithRT=General.concatenate(pair.getAuxScores(), filter.getAdditionalScores(pair));
				newResult.addStripe(pair.getPrimaryScore(), scoresWithRT, pair.getDeltaPrecursorMass(), pair.getDeltaFragmentMass(), pair.getMSMS());
				anyFoundWithRTFilter=true;
			}
		}
		
		// if nothing passes the RT filter then use the top match
		if (!anyFoundWithRTFilter) {
			float[] scoresWithRT=General.concatenate(startingBest.getAuxScores(), filter.getAdditionalScores(startingBest));
			newResult.addStripe(startingBest.getPrimaryScore(), scoresWithRT, startingBest.getDeltaPrecursorMass(), startingBest.getDeltaFragmentMass(), startingBest.getMSMS());
		}
		
		return newResult;
	}
	
	public int size() {
		return goodStripes.size();
	}
	
	/**
	 * assumes scores are already sorted
	 */
	public LibraryEntry getEntry() {
		/* assumes sorted in order, first is best scoring! */
		if (hasScoredResults()) {
			return goodStripes.get(0).getLibraryEntry();
		} else {
			return null;
		}
	}

	public void addStripe(float score, float[] auxScoreArray, float deltaPrecursorMass, float deltaFragmentMass, FragmentScan stripe) {
		throw new EncyclopediaException("TDC requires scoring library entries!");
	}

	public void addStripe(LibraryEntry entry, float score, float[] auxScoreArray, float deltaPrecursorMass, float deltaFragmentMass, FragmentScan stripe) {
		goodStripes.add(new ScoredPSM(entry, stripe, score, auxScoreArray, deltaPrecursorMass, deltaFragmentMass));
	}
	
	public void trim(int trimToN) {
		Collections.sort(goodStripes);
		Collections.reverse(goodStripes);
		while (goodStripes.size()>trimToN) {
			goodStripes.remove(goodStripes.size()-1);
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
		for (ScoredPSM pair : goodStripes) {
			scores.add(pair.getPrimaryScore());
		}
		float[] sorted=scores.toArray();
		Arrays.sort(sorted);
		return sorted;
	}
	
	public void setTrace(XYTraceInterface trace) {
		throw new EncyclopediaException("Trace not saved");
	}
	
	public XYTraceInterface getTrace() {
		throw new EncyclopediaException("Trace not recorded");
	}
	
	public boolean hasScoredResults() {
		return goodStripes.size()>0&&goodStripes.get(0).getMSMS()!=null;
	}
	
	public ScoredPSM getScoredMSMS() {
		if (hasScoredResults()) {
			return goodStripes.get(0);
		} else {
			return null;
		}
	}
	
	public ArrayList<ScoredPSM> getGoodMSMSCandidates() {
		return goodStripes;
	}
}

