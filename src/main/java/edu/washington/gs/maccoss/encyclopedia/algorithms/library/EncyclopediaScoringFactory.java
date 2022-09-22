package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;

public class EncyclopediaScoringFactory {
	public static final String V1_MODE_ARG = "-v1scoring";
	public static final String V2_MODE_ARG = "-v2scoring";
	public static boolean USE_LEGACY_SCORING_SYSTEM=false;
	
	public static LibraryScoringFactory getScoringFactory(HashMap<String, String> arguments, SearchParameters parameters) {
		LibraryScoringFactory factory;
		if (arguments.containsKey(V1_MODE_ARG)) {
			factory=new EncyclopediaOneScoringFactory(parameters);
		} else if (arguments.containsKey(V2_MODE_ARG)) {
			factory=new EncyclopediaTwoScoringFactory(parameters);
		} else {
			factory=getDefaultScoringFactory(parameters);
		}
		return factory;
	}
	public static LibraryScoringFactory getDefaultScoringFactory(SearchParameters parameters) {
		LibraryScoringFactory factory;
		if (USE_LEGACY_SCORING_SYSTEM) {
			factory=new EncyclopediaOneScoringFactory(parameters);
		} else {
			factory=new EncyclopediaTwoScoringFactory(parameters);
		}
		return factory;
	}
}
