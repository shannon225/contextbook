package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;

public class LibraryBackground {
	private final int[] background=new int[4000];
	private final int total;
	
	public static void main(String[] args) throws Exception {
		File libraryFile=new File("/Users/searleb/Documents/school/projects/pecandata/cptac2_human_hcd_selected.elib");
		LibraryFile file=new LibraryFile();
		file.openFile(libraryFile);
		ArrayList<LibraryEntry> entries=file.getEntries(new Range(700, 800));
		LibraryBackground background=new LibraryBackground(entries);
		for (int i=0; i<background.background.length; i++) {
			System.out.println(i+"\t"+background.getFraction(i+0.1)); // to avoid rounding errors
		}
	}

	public LibraryBackground(ArrayList<LibraryEntry> entries) {
		Arrays.fill(background, 1); // add a pseudocount
		for (LibraryEntry entry : entries) {
			double[] masses=entry.getMassArray();
			for (double mass : masses) {
				int index=(int)mass; // truncate
				if (index<background.length) {
					background[index]++;
				}
			}
		}
		int t=0;
		for (int i=0; i<background.length; i++) {
			t+=background[i];
		}
		this.total=t;
	}
	
	public float getFraction(double mass) {
		int index=(int)mass; // truncate
		int count=index>=background.length?1:background[index];
		return (total/(float)background.length)/(float)count;
	}
}
