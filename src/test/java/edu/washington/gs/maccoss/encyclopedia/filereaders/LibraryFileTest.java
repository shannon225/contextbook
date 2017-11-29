package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;

public class LibraryFileTest {
	public static void main(String[] args) throws Exception {
		int[] counts=new int[100];
		
		addCounts(counts, new File("/Volumes/BriansSSD/hela_serum_timecourse_library/23aug2017_hela_serum_timecourse_4mz_narrow_1.mzML.elib"));
		addCounts(counts, new File("/Volumes/BriansSSD/hela_serum_timecourse_library/23aug2017_hela_serum_timecourse_4mz_narrow_2.mzML.elib"));
		addCounts(counts, new File("/Volumes/BriansSSD/hela_serum_timecourse_library/23aug2017_hela_serum_timecourse_4mz_narrow_3.mzML.elib"));
		addCounts(counts, new File("/Volumes/BriansSSD/hela_serum_timecourse_library/23aug2017_hela_serum_timecourse_4mz_narrow_4.mzML.elib"));
		addCounts(counts, new File("/Volumes/BriansSSD/hela_serum_timecourse_library/23aug2017_hela_serum_timecourse_4mz_narrow_5.mzML.elib"));
		addCounts(counts, new File("/Volumes/BriansSSD/hela_serum_timecourse_library/23aug2017_hela_serum_timecourse_4mz_narrow_6.mzML.elib"));
		for (int i=0; i<counts.length; i++) {
			System.out.println(i+"\t"+counts[i]);
		}
		
		/*
		int[] counts=new int[100];
		
		addCounts(counts, new File("/Users/searleb/Documents/chromatogram_library_manuscript/quant_replicates/23aug2017_hela_serum_timecourse_pool_wide_001_170829031834.mzML.elib"));
		addCounts(counts, new File("/Users/searleb/Documents/chromatogram_library_manuscript/quant_replicates/23aug2017_hela_serum_timecourse_pool_wide_002.mzML.elib"));
		addCounts(counts, new File("/Users/searleb/Documents/chromatogram_library_manuscript/quant_replicates/23aug2017_hela_serum_timecourse_pool_wide_003.mzML.elib"));
		for (int i=0; i<counts.length; i++) {
			System.out.println(i+"\t"+counts[i]);
		}
		 */
	}

	public static void addCounts(int[] counts, File f) throws IOException, SQLException, DataFormatException {
		LibraryFile library=(LibraryFile)BlibToLibraryConverter.getFile(f);

		ArrayList<LibraryEntry> entries=library.getAllEntries(false);
		for (LibraryEntry entry : entries) {
			float[] c=entry.getCorrelationArray();
			int n=0;
			for (int i = 0; i < c.length; i++) {
				if (c[i]>=TransitionRefiner.quantitativeCorrelationThreshold) {
					n++;
				}
			}
			counts[n]++;
		}
	}
}
