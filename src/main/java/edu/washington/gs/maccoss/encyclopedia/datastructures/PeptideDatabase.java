package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.VariantFastaPeptideEntry;

public class PeptideDatabase implements Iterable<FastaPeptideEntry> {
	private final HashMap<String, FastaPeptideEntry> peptidesBySequence=new HashMap<String, FastaPeptideEntry>();
	
	public PeptideDatabase() {
	}
	
	/**
	 * note this is destructive to the FastaPeptideEntries (it can modify their annotation lists)
	 * @param newPeptide
	 */
	public void add(FastaPeptideEntry newPeptide) {
		FastaPeptideEntry entry=peptidesBySequence.get(newPeptide.getSequence());
		if (entry!=null) {
			if (entry instanceof VariantFastaPeptideEntry&&!(newPeptide instanceof VariantFastaPeptideEntry)) {
				// prefer to keep non-variant peptides as the canonical annotation
				peptidesBySequence.put(newPeptide.getSequence(), newPeptide);
				newPeptide.addAccessions(entry.getAccessions());
			} else {
				// if neither or both are variants then keep first
				entry.addAccessions(newPeptide.getAccessions());
			}
		} else {
			peptidesBySequence.put(newPeptide.getSequence(), newPeptide);
		}
	}
	
	public ArrayList<FastaPeptideEntry> getPeptides() {
		ArrayList<FastaPeptideEntry> entries=new ArrayList<FastaPeptideEntry>(peptidesBySequence.values());
		Collections.sort(entries);
		return entries;
	}
	
	@Override
	public Iterator<FastaPeptideEntry> iterator() {
		return peptidesBySequence.values().iterator();
	}
	
	public int size() {
		return peptidesBySequence.size();
	}
}
