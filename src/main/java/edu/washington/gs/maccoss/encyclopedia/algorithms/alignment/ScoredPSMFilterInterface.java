package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ScoredPSM;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;

public interface ScoredPSMFilterInterface {

	boolean passesFilter(ScoredPSM psm);

	float[] getAdditionalScores(ScoredPSM psm);

	void makePlots(SearchParameters params, ArrayList<ScoredPSM> psms, Optional<File> saveFileSeed);
}