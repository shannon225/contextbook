package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.List;

public interface ProteinGroupInterface extends Comparable<ProteinGroupInterface> {

	List<String> getSequences();

	List<String> getEquivalentAccessions();

	boolean isDecoy();
	
	float getNSPScore();
	
	float getPosteriorErrorProb();
}