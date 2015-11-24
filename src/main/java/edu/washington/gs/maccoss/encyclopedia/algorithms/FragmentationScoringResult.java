package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class FragmentationScoringResult extends PeptideScoringResult {
	ArrayList<XYTrace> fragmentationTraces=new ArrayList<XYTrace>();

	public void addFragmentationTrace(XYTrace trace) {
		fragmentationTraces.add(trace);
	}
	
	public ArrayList<XYTrace> getFragmentationTraces() {
		return fragmentationTraces;
	}
}
