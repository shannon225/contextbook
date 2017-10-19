package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;

public class ExtendedFastaEntry extends FastaEntry {
	private final ArrayList<AlleleVariant> potentialVariants=new ArrayList<>();

	public ExtendedFastaEntry(String filename, String annotation, String sequence) {
		super(filename, annotation, sequence);
	}

	public void addPotentialVariant(AlleleVariant variant) {
		potentialVariants.add(variant);
	}
}
