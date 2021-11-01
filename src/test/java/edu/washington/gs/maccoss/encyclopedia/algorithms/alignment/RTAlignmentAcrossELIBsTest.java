package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;
import gnu.trove.procedure.TObjectIntProcedure;

public class RTAlignmentAcrossELIBsTest {
	public static void main(String[] args) throws Exception {
		File[] searches=new File[] {
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2017apr30/20170430_HeLa_DIA_400_1000_60min_01.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2017apr30/20170430_HeLa_DIA_A_01_170506084818.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2017apr30/20170430_HeLa_pool_120min_4mz_DIA_frac_03.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2017apr30/20170430_HeLa_pool_120min_4mz_nonlinear_DIA_frac_03.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2017apr30/20170430_HeLa_DIA_B_01.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2017aug23/23aug2017_hela_serum_timecourse_pool_wide_003.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2017aug23/23aug2017_hela_serum_timecourse_pool_wide_002.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2017aug23/23aug2017_hela_serum_timecourse_pool_wide_001_170829031834.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2015dec11/121115_bcs_hela_24mz_400_1000_4c_0D_151220133553.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2015dec11/121115_bcs_hela_24mz_400_1000_m20c_0D.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2015dec11/121115_bcs_hela_24mz_400_1000_4c_0D.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2020jul29/2020jul29_BCS_explo_hela_15k_8mzol_90min_9to32.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2020jul29/2020jul29_BCS_explo_hela_15k_8mzol_60min.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2020jul29/2020jul29_BCS_explo_hela_15k_8mzol_120min.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2020jul29/2020jul29_BCS_explo_hela_15k_8mzol_90min_5to40.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2020jan16/bcs_2020jan16_hela_curve_A_9.dia.elib"),
				new File("/Users/searleb/Documents/OSU/damien/hela_raws/2020jan16/bcs_2020jan16_hela_curve_B_9.dia.elib")
		};

		File selection=new File("/Users/searleb/Documents/OSU/damien/hela_raws/detected_peptides.dlib");
		final TObjectFloatHashMap<String> rtsByPeptideModSeq = getRTs(selection);
		HashMap<File, TObjectFloatHashMap<String>> individualFiles=new HashMap<>();
		for (File file : searches) {
			individualFiles.put(file, getRTs(file));
		}
		TObjectIntHashMap<String> peptideFrequency=new TObjectIntHashMap<>();
		for (TObjectFloatHashMap<String> individual : individualFiles.values()) {
			for (String peptideModSeq : individual.keySet()) {
				peptideFrequency.adjustOrPutValue(peptideModSeq, 1, 1);
			}
		}
		
		HashSet<String> keepers=new HashSet<>();
		peptideFrequency.forEachEntry(new TObjectIntProcedure<String>() {
			@Override
			public boolean execute(String a, int b) {
				if (b>=10) {
					keepers.add(a);
				}
				return true;
			}
		});
		
		for (Entry<File, TObjectFloatHashMap<String>> entry : individualFiles.entrySet()) {
			final ArrayList<XYPoint> points=new ArrayList<>();

			entry.getValue().forEachEntry(new TObjectFloatProcedure<String>() {
				@Override
				public boolean execute(String a, float b) {
					float std=rtsByPeptideModSeq.get(a);
					points.add(new RTRTPoint(b/60f, std, false, a));
					return true;
				}
			});
			
			RetentionTimeFilter filter=RetentionTimeFilter.getFilter(points, "actual retention time (min)", "iRT");
			filter.plot(points, Optional.of(entry.getKey()), "actualRT", "iRT");
		}
		
		System.out.println(rtsByPeptideModSeq.size()+" vs "+keepers.size());
	}
	private static TObjectFloatHashMap<String> getRTs(File selection) throws IOException, SQLException, DataFormatException {
		LibraryFile select=new LibraryFile();
		select.openFile(selection);
		
		ArrayList<LibraryEntry> entries=select.getAllEntries(false, new AminoAcidConstants());
		TObjectFloatHashMap<String> rtsByPeptideModSeq=new TObjectFloatHashMap<>();
		for (LibraryEntry entry : entries) {
			rtsByPeptideModSeq.put(entry.getPeptideModSeq(), entry.getRetentionTime());
		}
		select.close();
		return rtsByPeptideModSeq;
	}
	public static void truncateLibraryToDetectedPeptides(String[] args) throws Exception {
		File primary=new File("/Users/searleb/Documents/OSU/damien/hela_raws/hela_specific.z2_nce33.dlib");
		LibraryFile lib=new LibraryFile();
		lib.openFile(primary);
		
		File selection=new File("/Users/searleb/Documents/OSU/damien/hela_raws/hela_specific_clib.z2_nce33.elib");
		LibraryFile select=new LibraryFile();
		select.openFile(selection);
		
		ArrayList<LibraryEntry> entries=select.getAllEntries(false, new AminoAcidConstants());
		HashSet<String> peptideModSeqs=new HashSet<>();
		for (LibraryEntry entry : entries) {
			peptideModSeqs.add(entry.getPeptideModSeq());
		}
		
		ArrayList<LibraryEntry> selectedEntries=new ArrayList<>();
		entries=lib.getAllEntries(false, new AminoAcidConstants());
		for (LibraryEntry entry : entries) {
			if (peptideModSeqs.contains(entry.getPeptideModSeq())) {
				selectedEntries.add(entry);
			}
		}

		LibraryFile result=new LibraryFile();
		result.openFile();
		result.dropIndices();
		result.addEntries(selectedEntries);
		result.addProteinsFromEntries(selectedEntries);
		result.createIndices();
		result.saveAsFile(new File("/Users/searleb/Documents/OSU/damien/hela_raws/detected_peptides.dlib"));
		result.close();
	}
}
