package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

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
			entry.addAccessions(newPeptide.getAccessions());
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
