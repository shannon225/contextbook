package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class IARPATestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false, true, false);
	
	private static final File[] files=new File[] {
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_11_1.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_12_2.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_13_3.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_16_4.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_17_5.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_18_6.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_19_7.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_20_8.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_23_9.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_24_12.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_25_13.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_26_14.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_27_15.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_30_16.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_31_17.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_32_18.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_33_19.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_34_20.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_46_21.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_47_22.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_48_23.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_49_24.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_50_25.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_53_26.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_54_27.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_55_28.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_56_29.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_57_30.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_60_31.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_61_32.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_62_33.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_63_34.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_64_35.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_67_36.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_68_37.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_69_38.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_70_39.dia.elib"),
			new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/XXX_2019_0304_RJ_71_40.dia.elib")
	};

	public static void main(String[] args) throws Exception {
		File userFile = new File("/Users/searleb/Documents/iarpa/individuals/localized_individual_results/library/clib.elib");
		TObjectFloatHashMap<String> narrow = getRTs(userFile);
		
		HashSet<String> allPeptides=new HashSet<>();
		HashMap<File, TObjectFloatHashMap<String>> alignments=new HashMap<>();
		for (File f : files) {
			TObjectFloatHashMap<String> alignment = getAlignment(f, narrow);
			alignments.put(f, alignment);
			allPeptides.addAll(alignment.keySet());
		}

		System.out.print("peptide");
		for (File f : files) {
			System.out.print("\t"+f.getName());			
		}
		System.out.println();
		
		for (String peptide : allPeptides) {
			System.out.print(peptide);
			for (File f : files) {
				System.out.print("\t");
				
				float alignedRT=alignments.get(f).get(peptide);
				System.out.print(alignedRT);
			}
			System.out.println();
		}
	}

	private static TObjectFloatHashMap<String> getAlignment(File f, TObjectFloatHashMap<String> narrow)
			throws IOException, SQLException, DataFormatException {
		TObjectFloatHashMap<String> wide = getRTs(f);
		
		ArrayList<XYPoint> points=new ArrayList<XYPoint>();
		
		wide.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String a, float b) {
				float alt=narrow.get(a);
				if (narrow.getNoEntryValue()!=alt) {
					// narrow first, then wide
					points.add(new XYPoint(alt, b));
				}
				return true;
			}
		});
		
		RetentionTimeAlignmentInterface alignment=RetentionTimeFilter.getFilter(points, "narrow", f.getName(), 10000);
		
		TObjectFloatHashMap<String> aligned = new TObjectFloatHashMap<>(); 
		wide.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String a, float b) {
				float alt=alignment.getXValue(b);
				aligned.put(a, alt);
				return true;
			}
		});
		return aligned;
	}

	private static TObjectFloatHashMap<String> getRTs(File userFile)
			throws IOException, SQLException, DataFormatException {
		LibraryFile lf=new LibraryFile();
		lf.openFile(userFile);
		
		ArrayList<LibraryEntry> entries=lf.getAllEntries(false, PARAMETERS.getAAConstants());
		TObjectFloatHashMap<String> rtInSec=new TObjectFloatHashMap<>();
		for (LibraryEntry entry : entries) {
			rtInSec.put(entry.getPeptideModSeq(), entry.getRetentionTime());
		}
		return rtInSec;
	}
}
