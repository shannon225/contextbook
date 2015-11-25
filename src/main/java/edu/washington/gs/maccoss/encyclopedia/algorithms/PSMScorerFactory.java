package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;

public interface PSMScorerFactory {
	public PSMScorer getScorer(SearchParameters parameters);
}
