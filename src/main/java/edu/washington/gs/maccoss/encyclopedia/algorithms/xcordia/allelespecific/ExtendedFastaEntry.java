package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;

public class ExtendedFastaEntry extends FastaEntry {
	ArrayList<VariableModification> potentialMods=new ArrayList<>();

	public ExtendedFastaEntry(String filename, String annotation, String sequence) {
		super(filename, annotation, sequence);
	}

}
