package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Optional;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeptideXYPoint;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RTAlignmentAcrossELIBsTest;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class LibraryUtilitiesTest {
	public static void main(String[] args) throws Exception {
		File dir=new File("/Volumes/MacOnlySSD/2023_06_29_CD8_TCell_Exhaustion_Data/dia_files/libraries/");
		File lib1=new File(dir, "051622_Mouse_Tcell_pool_clib.elib");
		File lib2=new File(dir, "combined_2x_pooled_GPF_day8_tcells.dia.elib");
		File lib3=new File(dir, "lcmv_6x_GPF_library_400-1000.elib");
		File lib4=new File(dir, "global_invivo_invitro_pooled_library.elib");
		
		String[] targetAccessions = new String[] { "E9PVX6", "O35284", "O54839", "O54890", "O70309", "O88508", "O88713",
				"O88875", "P01580", "P01590", "P04187", "P05555", "P06876", "P07901", "P08113", "P09055", "P09793",
				"P10820", "P11499", "P11688", "P11835", "P13864", "P14211", "P14602", "P15379", "P16045", "P16627",
				"P16872", "P17156", "P17879", "P18337", "P19324", "P20029", "P23611", "P24063", "P26011", "P27782",
				"P28033", "P31041", "P32972", "P35385", "P37217", "P38647", "P41183", "P43406", "P48722", "P50283",
				"P55772", "P56528", "P58252", "P58681", "P60843", "P63017", "P63037", "P63038", "P63073", "P86176",
				"P97303", "Q00417", "Q00651", "Q02242", "Q03347", "Q04683", "Q07763", "Q3TTA7", "Q5EBG6", "Q60636",
				"Q60765", "Q61024", "Q61316", "Q61696", "Q61699", "Q61735", "Q61739", "Q61790", "Q64131", "Q64433",
				"Q66JW3", "Q6NZJ6", "Q6ZWX6", "Q8BK58", "Q8BM72", "Q8K0U4", "Q8R180", "Q8VIM0", "Q91VC3", "Q91YW3",
				"Q921X9", "Q922R8", "Q99M31", "Q99MB1", "Q99P31", "Q99PR8", "Q9CQN1", "Q9CZJ2", "Q9D659", "Q9D6H2",
				"Q9DAM3", "Q9EQY0", "Q9ET39", "Q9JJK5", "Q9JK92", "Q9JKD8", "Q9JKR6", "Q9JLV1", "Q9QUN7", "Q9QXT0",
				"Q9QYJ3", "Q9QZ57", "Q9R0G7", "Q9R1E0", "Q9WVS0", "Q9Z0T9", "Q9Z2B5", };

		AminoAcidConstants aaConstants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());

		TObjectFloatHashMap<String> rtsByPeptideModSeq=RTAlignmentAcrossELIBsTest.getRTs(lib1);
		HashMap<String, LibraryEntry> entriesByPeptideModSeq=new HashMap<String, LibraryEntry>();
		HashMap<File, RetentionTimeFilter> filtersByFile=new HashMap<>();

		LibraryFile library=new LibraryFile();
		library.openFile(lib1);
		ArrayList<LibraryEntry> entries=library.getAllEntries(false, aaConstants);
		EACHENTRY: for (LibraryEntry entry : entries) {
			for (String accession : entry.getAccessions()) {
				for (int i = 0; i < targetAccessions.length; i++) {
					if (accession.indexOf(targetAccessions[i])>=0) {
						entriesByPeptideModSeq.put(entry.getPeptideModSeq(), entry);
						continue EACHENTRY;
					}
				}
			}
		}
		library.close();
		System.out.println(lib1.getName()+" --> "+entriesByPeptideModSeq.size());
		
		for (File f : new File[] {lib2, lib3, lib4}) {
			TObjectFloatHashMap<String> altRTs=RTAlignmentAcrossELIBsTest.getRTs(f);

			final ArrayList<XYPoint> points=new ArrayList<>();

			altRTs.forEachEntry(new TObjectFloatProcedure<String>() {
				@Override
				public boolean execute(String a, float b) {
					float std=rtsByPeptideModSeq.get(a);
					if (std>0) {
						points.add(new PeptideXYPoint(b/60f, std/60f, false, a));
					}
					return true;
				}
			});
			
			RetentionTimeFilter filter=RetentionTimeFilter.getFilter(points, f.getName(), lib1.getName());
			filter.plot(points, Optional.empty(), f.getName(), lib1.getName());
			filtersByFile.put(f, filter);

			library=new LibraryFile();
			library.openFile(f);
			entries=library.getAllEntries(false, aaConstants);
			EACHENTRY: for (LibraryEntry entry : entries) {
				for (String accession : entry.getAccessions()) {
					for (int i = 0; i < targetAccessions.length; i++) {
						if (accession.indexOf(targetAccessions[i])>=0) {
							LibraryEntry updated=entry.updateRetentionTime(60f*filter.getYValue(entry.getRetentionTime()/60f));
							entriesByPeptideModSeq.put(entry.getPeptideModSeq(), updated);
							continue EACHENTRY;
						}
					}
				}
			}
			library.close();
			System.out.println(f.getName()+" --> "+entriesByPeptideModSeq.size());
		}

		ArrayList<LibraryEntry> targets=new ArrayList<>(entriesByPeptideModSeq.values());
		
		for (File f : new File[] {lib1, lib2, lib3, lib4}) {
			File libFile=new File(dir, "targets_"+f.getName());
			library=new LibraryFile();
			library.openFile();
			library.dropIndices();

			if (filtersByFile.containsKey(f)) {
				RetentionTimeFilter filter=filtersByFile.get(f);
				ArrayList<LibraryEntry> updatedTargets=new ArrayList<>();
				for (LibraryEntry entry : targets) {
					LibraryEntry updated=entry.updateRetentionTime(60f*filter.getXValue(entry.getRetentionTime()/60f));
					updatedTargets.add(updated);
				}
				library.addEntries(updatedTargets);
				
			} else {
				library.addEntries(targets);
			}
			
			library.addProteinsFromEntries(entries);
			library.createIndices();
			library.saveAsFile(libFile);
			library.close();
		}
		
	}
	
	public static void main4(String[] args) throws Exception {
		File inFile=new File("/Volumes/bcsbluessd/TPAD/combined_library_noquant.elib");
		File saveFile=new File("/Volumes/bcsbluessd/TPAD/Vpool_combined_library_noquant.dlib");
		File rtFile=new File("/Volumes/bcsbluessd/TPAD/Vpool_combined_library_quant.dlib");

		// put all RTs into map from Vpool_combined_library_quant.dlib (the extracted RTs from the global quant analysis)
		AminoAcidConstants aaConstants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
		LibraryFile library=new LibraryFile();
		library.openFile(rtFile);
		TObjectFloatHashMap<String> rtMap=new TObjectFloatHashMap<>();
		for (LibraryEntry entry : library.getAllEntries(false, aaConstants)) {
			String key=entry.getPeptideModSeq()+"+"+entry.getPrecursorCharge();
			rtMap.put(key, entry.getScanStartTime());
		}
		library.close();
		
		library=new LibraryFile();
		library.openFile(inFile);

		// for each entry in /Volumes/bcsbluessd/TPAD/combined_library_noquant.elib
		HashMap<String, LibraryEntry> toWrite=new HashMap<>();
		for (LibraryEntry entry : library.getAllEntries(false, aaConstants)) {
			String key=entry.getPeptideModSeq()+"+"+entry.getPrecursorCharge();
			if (rtMap.contains(key)) {
				// if we have an RT form the global quant analysis, then update the RT if the combined quant RT is better
				if (toWrite.containsKey(key)) {
					if (toWrite.get(key).getScore()>entry.getScore()) {
						toWrite.put(key, entry.updateRetentionTime(rtMap.get(key)));
					}
				} else {
					toWrite.put(key, entry.updateRetentionTime(rtMap.get(key)));
				}
			}
		}
		ArrayList<LibraryEntry> entries=new ArrayList<>(toWrite.values());
		Collections.sort(entries);
		
		Logger.logLine("Found "+toWrite.size()+" peptides. Writing to ["+saveFile.getAbsolutePath()+"]...");

		LibraryFile saveLibrary=new LibraryFile();
		saveLibrary.openFile();
		saveLibrary.dropIndices();
		saveLibrary.addEntries(entries);
		saveLibrary.addProteinsFromEntries(entries);
		saveLibrary.createIndices();
		saveLibrary.saveAsFile(saveFile);
		
		saveLibrary.close();
	}
	
	public static void main7(String[] args) throws Exception {
		File lib=new File("/Users/searleb/Downloads/lit_dilution_ms2.dlib");

		LibraryFile library=new LibraryFile();
		library.openFile(lib);

		ArrayList<String> sourceFiles=new ArrayList<String>();
		Connection c=library.getConnection();
		try {
			Statement s=c.createStatement();
			try {
				
				Logger.logLine("Getting source files...");
				ResultSet rs=s.executeQuery("select distinct SourceFile from entries");
				while (rs.next()) {
					sourceFiles.add(rs.getString(1));
				}
				rs.close();
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
		
		for (String source : sourceFiles) {
			String sub=source.substring("20221222_HB_TP_Evo_Whisper100_40SPD_Hefe_".length(), source.length()-".raw".length());
			System.out.println(source+"\t"+sub);
			File saveFile=new File(lib.getParentFile(), sub+".dlib");
			HashSet<String> targets=new HashSet<String>();
			targets.add(source);
			LibraryUtilities.subsetLibrary(saveFile, 0, 1000*60, 0, 10000, targets, library);
		}

		library.close();
	}
	
	public static void main2(String[] args) throws Exception {

		File twoDLC=new File("/Users/searleb/Downloads/msms.dlib");

		LibraryFile library=new LibraryFile();
		library.openFile(twoDLC);

		float[] timeBoundaries=new float[] {0f, 2f, 32f, 62f, 92f, 122f, 152f, 182f, 212f, 242f, 272f, 302f, 332f, 362f, 392f, 1000f};
		for (int i=1; i<timeBoundaries.length; i++) {
			File saveFile=new File(twoDLC.getParentFile(), "msms_"+Math.round(timeBoundaries[i-1])+"to"+Math.round(timeBoundaries[i])+".dlib");
			System.out.println("writing "+saveFile.getName());
			LibraryUtilities.subsetLibrary(saveFile, timeBoundaries[i-1]*60f, timeBoundaries[i]*60f, 0, 100000, new HashSet<>(), library);
		}

		library.close();
	}
	
	public static void main3(String[] args) throws Exception {
		File inFile = new File("/Users/searleb/Documents/cobbs/2020dec30_cobbs/2020dec03_cobbs_cmv_inf_gpfdia_clib.elib");
		File saveFile = new File("/Users/searleb/Documents/cobbs/2020dec30_cobbs/2020dec03_cobbs_cmv_inf_gpfdia_clib_hcmv_only.elib");
		LibraryFile library=new LibraryFile();
		library.openFile(inFile);
		String accessionNumberKeyword="HCMV";

		LibraryFile saveLibrary=new LibraryFile();
		saveLibrary.openFile();
		
		ArrayList<LibraryEntry> toWrite=new ArrayList<>();
		for (LibraryEntry entry : library.getAllEntries(false, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()))) {
			for (String accession : entry.getAccessions()) {
				if (accession.contains(accessionNumberKeyword)) {
					toWrite.add(entry);
				}
			}
		}
		Logger.logLine("Found "+toWrite.size()+" peptides for "+accessionNumberKeyword+". Writing to ["+saveFile.getAbsolutePath()+"]...");
		
		saveLibrary.dropIndices();
		saveLibrary.addEntries(toWrite);
		saveLibrary.addProteinsFromEntries(toWrite);
		saveLibrary.createIndices();
		saveLibrary.saveAsFile(saveFile);
		
		saveLibrary.close();
	}
	
	public static void main6(String[] args) throws Exception {
		//File dirFile=new File("/Users/searleb/Documents/damien/dda_library_search/non-tryptics/");
		File dirFile=new File("/Users/searle.30/Downloads/");
		File[] inFiles=new File[] {
				//new File(dirFile, "uniprot_human-reference_reviewed_2022mar02.prosit_input.aspn_nce29.prosit_cid2020.dlib"),
				//new File(dirFile, "uniprot_human-reference_reviewed_2022mar02.prosit_input.aspn_nce29.prosit_hcd2020.dlib"),
				//new File(dirFile, "uniprot_human-reference_reviewed_2022mar02.prosit_input.gluc_nce29.prosit_cid2020.dlib"),
				//new File(dirFile, "uniprot_human-reference_reviewed_2022mar02.prosit_input.gluc_nce29.prosit_hcd2020.dlib"),
				new File(dirFile, "uniprot_human-reference_reviewed_2022mar02.prosit_input.tryp_nce29.prosit_cid2020.dlib"),
				new File(dirFile, "uniprot_human-reference_reviewed_2022mar02.prosit_input.tryp_nce29.prosit_hcd2020.dlib"),
				new File(dirFile, "uniprot_human-reference_reviewed_2022mar02.prosit_input.tryp_nce34.prosit_hcd2020.dlib")};
		for (File inFile : inFiles) {
			File saveFile=new File(dirFile, "mzcorrected_"+inFile.getName());
			//File inFile=new File("/Users/searleb/Documents/damien/dda_library_search/non-tryptics/scribe_regression/uniprot_human-reference_reviewed_2022mar02.prosit_input.tryp_nce29.prosit_hcd2020.dlib");	
			//File saveFile=new File("/Users/searleb/Documents/damien/dda_library_search/non-tryptics/scribe_regression/mzcorrected_uniprot_human-reference_reviewed_2022mar02.prosit_input.tryp_nce29.prosit_hcd2020.dlib");
			AminoAcidConstants constants=new AminoAcidConstants();
			
			LibraryFile library=new LibraryFile();
			library.openFile(inFile);
	
			LibraryFile saveLibrary=new LibraryFile();
			saveLibrary.openFile();
			saveLibrary.dropIndices();
			
			HashMap<String, HashSet<String>> targetAccessionsByPeptide=new HashMap<>();
			HashMap<String, HashSet<String>> decoyAccessionsByPeptide=new HashMap<>();
			
			ArrayList<LibraryEntry> toWrite=new ArrayList<>();
			TDoubleArrayList deltamz=new TDoubleArrayList();
			for (int i = 0; i <= 10; i++) {
				float min=i*100;
				float max=i==10?10000.0f:((i+1)*100);
				Range range=new Range(min, max);
				System.out.print("batch "+range.toString());
			
				for (LibraryEntry entry : library.getEntries(range, false, constants)) {
					double mz=constants.getChargedMass(entry.getPeptideModSeq(), entry.getPrecursorCharge());
					deltamz.add(mz-entry.getPrecursorMZ());
					toWrite.add(entry.updatePrecursorMz(mz));
	
					HashMap<String, HashSet<String>> map;
					if (entry.isDecoy()) {
						map=decoyAccessionsByPeptide;
					} else {
						map=targetAccessionsByPeptide;
					}
					HashSet<String> accessions=map.get(entry.getPeptideSeq());
					if (accessions==null) {
						accessions=new HashSet<>();
						map.put(entry.getPeptideSeq(), accessions);
					}
					accessions.addAll(entry.getAccessions());
					
					if (toWrite.size()>10000) {
						saveLibrary.addEntries(toWrite);
						toWrite.clear();
						System.out.print('.');
					}
					
				}
	
				if (toWrite.size()>0) {
					saveLibrary.addEntries(toWrite);
					toWrite.clear();
					System.out.println("Finished batch.");
				}
			}
			Logger.logLine("Found "+toWrite.size()+" peptides. Writing to ["+saveFile.getAbsolutePath()+"]...");
			ArrayList<XYPoint> points=PivotTableGenerator.createPivotTable(General.toFloatArray(deltamz.toArray()));
			
			XYTrace xyTrace = new XYTrace(points, GraphType.line, "delta mz");
			System.out.println("Most common mass error: "+xyTrace.getMaxXY().x+" m/z");
			//Charter.launchChart("delta mz", "count", true, xyTrace);
			
			saveLibrary.addProteinsFromEntries(targetAccessionsByPeptide, decoyAccessionsByPeptide);
			saveLibrary.createIndices();
			saveLibrary.saveAsFile(saveFile);
			
			saveLibrary.close();
		}
	}
	
	public static void main5(String[] args) throws Exception {
		File inFile=new File("/Volumes/bcsbluessd/TPAD/combined_library_quant.elib");
		File saveFile=new File("/Volumes/bcsbluessd/TPAD/libs/");
		LibraryFile library=new LibraryFile();
		library.openFile(inFile);
		LibraryUtilities.extractSampleSpecificLibraries(new EmptyProgressIndicator(), saveFile, Optional.ofNullable((HashMap<String, double[]>)null), library, SearchParameterParser.getDefaultParametersObject());
		library.close();
	}
}
