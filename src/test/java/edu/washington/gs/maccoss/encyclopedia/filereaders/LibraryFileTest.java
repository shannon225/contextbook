package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class LibraryFileTest {
	public static void main(String[] args) throws Exception {
		LibraryInterface ddaLib=BlibToLibraryConverter.getFile(new File("/Users/searleb/Documents/projects/encyclopedia/HumanTotalProteome/libs/HeLa.dlib"));
		final ArrayList<LibraryEntry> ddaEntries=ddaLib.getAllEntries(false, new AminoAcidConstants());
		System.out.println("DDA: "+ddaEntries.size());
		
		LibraryInterface diaLib=BlibToLibraryConverter.getFile(new File("/Users/searleb/Documents/projects/encyclopedia/HumanTotalProteome/libs/22oct2017_hela_serum_timecourse_narrow_library.elib"));
		final ArrayList<LibraryEntry> diaEntries=diaLib.getAllEntries(false, new AminoAcidConstants());
		System.out.println("DIA: "+diaEntries.size());
		
		HashMap<String, LibraryEntry> ddaEntriesByPeptideModSeq=new HashMap<>();
		for (LibraryEntry entry : ddaEntries) {
			ddaEntriesByPeptideModSeq.put(getKey(entry), entry);
		}
		
		ArrayList<LibraryEntry> diaEntriesWithDDATimes=new ArrayList<>();
		ArrayList<LibraryEntry> diaEntriesWithDDAPeaks=new ArrayList<>();
		ArrayList<LibraryEntry> diaEntriesWithDDATimesAndPeaks=new ArrayList<>();
		
		for (LibraryEntry diaEntry : diaEntries) {
			LibraryEntry ddaEntry=ddaEntriesByPeptideModSeq.get(getKey(diaEntry));
			if (ddaEntry!=null) {
				diaEntriesWithDDATimes.add(diaEntry.updateRetentionTime(ddaEntry.getRetentionTime()));
				diaEntriesWithDDAPeaks.add(diaEntry.updateMS2(ddaEntry.getMassArray(), ddaEntry.getIntensityArray()));
				diaEntriesWithDDATimesAndPeaks.add(diaEntry.updateMS2(ddaEntry.getMassArray(), ddaEntry.getIntensityArray()).updateRetentionTime(ddaEntry.getRetentionTime()));
			}
		}
		
		LibraryFile diaLibraryWithDDARTsAndPeaks=new LibraryFile();
		diaLibraryWithDDARTsAndPeaks.openFile();
		diaLibraryWithDDARTsAndPeaks.addEntries(diaEntriesWithDDATimesAndPeaks);
		diaLibraryWithDDARTsAndPeaks.saveAsFile(new File("/Users/searleb/Documents/projects/encyclopedia/HumanTotalProteome/libs/chromatogramLibraryWithDDATimesAndPeaks.dlib"));
		diaLibraryWithDDARTsAndPeaks.close();

		LibraryFile diaLibraryWithDDATimes=new LibraryFile();
		diaLibraryWithDDATimes.openFile();
		diaLibraryWithDDATimes.addEntries(diaEntriesWithDDATimes);
		diaLibraryWithDDATimes.saveAsFile(new File("/Users/searleb/Documents/projects/encyclopedia/HumanTotalProteome/libs/chromatogramLibraryWithDDATimes.dlib"));
		diaLibraryWithDDATimes.close();
		
		LibraryFile diaLibraryWithDDAPeaks=new LibraryFile();
		diaLibraryWithDDAPeaks.openFile();
		diaLibraryWithDDAPeaks.addEntries(diaEntriesWithDDAPeaks);
		diaLibraryWithDDAPeaks.saveAsFile(new File("/Users/searleb/Documents/projects/encyclopedia/HumanTotalProteome/libs/chromatogramLibraryWithDDAPeaks.dlib"));
		diaLibraryWithDDAPeaks.close();
	}

	public static String getKey(LibraryEntry entry) {
		return entry.getPeptideModSeq()+"+"+entry.getPrecursorCharge();
	}
	
	public static void main2(String[] args) throws Exception {
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

		final AminoAcidConstants aaConstants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());

		ArrayList<LibraryEntry> entries=library.getAllEntries(false, aaConstants);
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
