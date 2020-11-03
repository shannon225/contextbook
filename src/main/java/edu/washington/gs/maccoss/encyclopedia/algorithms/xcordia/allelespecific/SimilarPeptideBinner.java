package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.util.ArrayList;
import java.util.Collection;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;

public class SimilarPeptideBinner {
	public SimilarPeptideBinner() {
	}
	
	public ArrayList<ArrayList<LibraryEntry>> binEntries(Collection<LibraryEntry> peptides) {
		ArrayList<ArrayList<LibraryEntry>> sets=new ArrayList<>();
		for (LibraryEntry peptide : peptides) {
			boolean added=false;
			for (ArrayList<LibraryEntry> arrayList : sets) {
				for (LibraryEntry fastaPeptideEntry : arrayList) {
					if (areSimilarEnough(peptide.getPeptideSeq(), fastaPeptideEntry.getPeptideSeq())) {
						added=true;
						break;
					}
				}
				if (added) {
					arrayList.add(peptide);
					break;
				}
			}
			if (!added) {
				ArrayList<LibraryEntry> newSet=new ArrayList<>();
				newSet.add(peptide);
				sets.add(newSet);
			}
		}
		return sets;
	}
	
	public ArrayList<ArrayList<FastaPeptideEntry>> binPeptides(Collection<FastaPeptideEntry> peptides) {
		ArrayList<ArrayList<FastaPeptideEntry>> sets=new ArrayList<>();
		for (FastaPeptideEntry peptide : peptides) {
			boolean added=false;
			for (ArrayList<FastaPeptideEntry> arrayList : sets) {
				for (FastaPeptideEntry fastaPeptideEntry : arrayList) {
					if (areSimilarEnough(peptide.getSequenceWithModsStripped(), fastaPeptideEntry.getSequenceWithModsStripped())) {
						added=true;
						break;
					}
				}
				if (added) {
					arrayList.add(peptide);
					break;
				}
			}
			if (!added) {
				ArrayList<FastaPeptideEntry> newSet=new ArrayList<>();
				newSet.add(peptide);
				sets.add(newSet);
			}
		}
		return sets;
	}
	
	private static final byte one=1;
	private static final byte two=2;
	public static boolean areSimilarEnough(String p1, String p2) {
		if (p1.length()==p2.length()) {
			byte deviantBits=0;
			for (int i = 0; i < p1.length(); i++) {
				if (!areSimilar(p1.charAt(i), p2.charAt(i))) {
					deviantBits+=one;
					if (deviantBits>two) return false; 
				}
			}
			return true;
		}
		return false;
	}
	
	public static boolean areSimilar(char a, char b) {
		if (a==b) return true;
		if ((a=='I'||a=='L')&&(b=='I'||b=='L')) {
			return true;
		}

		if ((a=='D'||a=='N')&&(b=='D'||b=='N')) {
			return true;
		}

		if ((a=='Q'||a=='E'||a=='K')&&(b=='Q'||b=='E'||b=='K')) {
			return true;
		}
		return false;
	}
}