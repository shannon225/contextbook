package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class LibraryFileTest {

	public static void main(String[] args) throws Exception {
		LibraryInterface lib=BlibToLibraryConverter.getFile(new File("/Volumes/bcsbluessd/kkolotyuk/inputs", "prosit-output-background.dlib"));
		final ArrayList<LibraryEntry> entries=lib.getEntries(new Range(339, 400), false, new AminoAcidConstants());
		System.out.println(entries.size());
		ArrayList<LibraryEntry> selected=new ArrayList<>();
		for (LibraryEntry entry : entries) {
			//if (entry.getPeptideSeq().startsWith("A")) {
				selected.add(entry);
			//}
		}
		System.out.println(selected.size());
		
		LibraryFile out=new LibraryFile();
		out.openFile();
		out.dropIndices();
		out.addEntries(selected);
		out.createIndices();
		out.saveAsFile(new File("/Volumes/bcsbluessd/kkolotyuk/inputs", "selected-prosit-output-background.dlib"));
		
	}
	
	public static void main3(String[] args) throws Exception {
		LibraryInterface ddaLib=BlibToLibraryConverter.getFile(new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce33.dlib"));
		final ArrayList<LibraryEntry> ddaEntries=ddaLib.getAllEntries(false, new AminoAcidConstants());
		System.out.println("DDA: "+ddaEntries.size());
		
		LibraryInterface diaLib=BlibToLibraryConverter.getFile(new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce33_clib.elib"));
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
		diaLibraryWithDDARTsAndPeaks.saveAsFile(new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce33_clibWithDDATimesAndPeaks.dlib"));
		diaLibraryWithDDARTsAndPeaks.close();

		LibraryFile diaLibraryWithDDATimes=new LibraryFile();
		diaLibraryWithDDATimes.openFile();
		diaLibraryWithDDATimes.addEntries(diaEntriesWithDDATimes);
		diaLibraryWithDDATimes.saveAsFile(new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce33_clibWithDDATimes.dlib"));
		diaLibraryWithDDATimes.close();
		
		LibraryFile diaLibraryWithDDAPeaks=new LibraryFile();
		diaLibraryWithDDAPeaks.openFile();
		diaLibraryWithDDAPeaks.addEntries(diaEntriesWithDDAPeaks);
		diaLibraryWithDDAPeaks.saveAsFile(new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce33_clibWithDDAPeaks.dlib"));
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
