package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.ScoredPSMFilterInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;

public abstract class AbstractScoringResult {
	public static final AbstractScoringResult POISON_RESULT=new PeptideScoringResult(null);
	
	
	public abstract AbstractScoringResult rescore(ScoredPSMFilterInterface filter);
	
	public abstract int size();
	
	public abstract LibraryEntry getEntry();

	public abstract void addStripe(float score, float[] auxScoreArray, float deltaPrecursorMass, float deltaFragmentMass, FragmentScan stripe);
	
	public abstract float getBestScore();
	public abstract float getSecondBestScore();
	
	public abstract void sort(int trimToN);
	
	public abstract void setTrace(XYTraceInterface trace);
	
	public abstract XYTraceInterface getTrace();
	
	public abstract boolean hasScoredResults();
	
	public abstract ScoredPSM getScoredMSMS();
	
	public abstract ArrayList<ScoredPSM> getGoodMSMSCandidates();
}
