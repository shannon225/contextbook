package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class FragmentationScoringResult extends PeptideScoringResult {
	ArrayList<XYTrace> fragmentationTraces=new ArrayList<XYTrace>();

	public FragmentationScoringResult(LibraryEntry entry) {
		super(entry);
	}

	public void addFragmentationTrace(XYTrace trace) {
		fragmentationTraces.add(trace);
	}

	public ArrayList<XYTrace> getFragmentationTraces() {
		return fragmentationTraces;
	}
}
