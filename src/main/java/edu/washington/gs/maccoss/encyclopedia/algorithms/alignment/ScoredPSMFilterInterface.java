package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;

public interface ScoredPSMFilterInterface {

	boolean passesFilter(PSMInterface psm);
	
	float getYRT(float xrt);

	float[] getAdditionalScores(PSMInterface psm);

	void makePlots(SearchParameters params, ArrayList<PSMInterface> psms, Optional<File> saveFileSeed);
}