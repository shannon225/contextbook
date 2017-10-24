package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;

public class LibraryFileTest {
	public static void main(String[] args) throws Exception {
		File f=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/hela/20170430_HeLa_phosp_DIA_B_01_170506220515.dia.elib");
		LibraryFile library=(LibraryFile)BlibToLibraryConverter.getFile(f);

		ArrayList<LibraryEntry> entries=library.getAllEntries(false);
		int count=0;
		for (LibraryEntry entry : entries) {
			float[] c=entry.getCorrelationArray();
			int n=0;
			for (int i = 0; i < c.length; i++) {
				if (c[i]>=TransitionRefiner.quantitativeCorrelationThreshold) {
					n++;
				}
			}
			if (n>=3) {
				count++;
			}
		}
		System.out.println(count+"/"+entries.size());
	}
}
