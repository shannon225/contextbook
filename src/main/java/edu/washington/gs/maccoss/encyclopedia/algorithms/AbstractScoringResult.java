package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.AbstractRetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TFloatArrayList;

public abstract class AbstractScoringResult {
	public static final AbstractScoringResult POISON_RESULT=new PeptideScoringResult(null);
	
	
	public abstract AbstractScoringResult rescore(RetentionTimeAlignmentInterface filter);
	
	public abstract int size();
	
	public abstract LibraryEntry getEntry();

	public abstract void addStripe(float score, float[] auxScoreArray, FragmentScan stripe);
	
	public abstract float getBestScore();
	public abstract float getSecondBestScore();
	
	public abstract void setTrace(XYTraceInterface trace);
	
	public abstract XYTraceInterface getTrace();
	
	public abstract boolean hasScoredResults();
	
	public abstract Pair<ScoredObject<FragmentScan>, float[]> getScoredMSMS();
	
	public abstract ArrayList<Pair<ScoredObject<FragmentScan>, float[]>> getGoodMSMSCandidates();
}
