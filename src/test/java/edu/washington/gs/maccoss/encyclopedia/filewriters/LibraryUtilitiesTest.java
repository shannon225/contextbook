package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class LibraryUtilitiesTest {
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
	
	public static void main(String[] args) throws Exception {
		File dirFile=new File("/Users/searleb/Documents/damien/dda_library_search/non-tryptics/");
		File[] inFiles=new File[] {new File(dirFile, "uniprot_human-reference_reviewed_2022mar02.prosit_input.aspn_nce29.prosit_cid2020.dlib"),
				new File(dirFile, "uniprot_human-reference_reviewed_2022mar02.prosit_input.aspn_nce29.prosit_hcd2020.dlib"),
				new File(dirFile, "uniprot_human-reference_reviewed_2022mar02.prosit_input.gluc_nce29.prosit_cid2020.dlib"),
				new File(dirFile, "uniprot_human-reference_reviewed_2022mar02.prosit_input.gluc_nce29.prosit_hcd2020.dlib"),
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
			
			ArrayList<LibraryEntry> toWrite=new ArrayList<>();
			TDoubleArrayList deltamz=new TDoubleArrayList();
			for (LibraryEntry entry : library.getAllEntries(false, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()))) {
				double mz=constants.getChargedMass(entry.getPeptideModSeq(), entry.getPrecursorCharge());
				deltamz.add(mz-entry.getPrecursorMZ());
				toWrite.add(entry.updatePrecursorMz(mz));
				
			}
			Logger.logLine("Found "+toWrite.size()+" peptides. Writing to ["+saveFile.getAbsolutePath()+"]...");
			ArrayList<XYPoint> points=PivotTableGenerator.createPivotTable(General.toFloatArray(deltamz.toArray()));
			//Charter.launchChart("delta mz", "count", true, new XYTrace(points, GraphType.line, "delta mz"));
			
			saveLibrary.dropIndices();
			saveLibrary.addEntries(toWrite);
			saveLibrary.addProteinsFromEntries(toWrite);
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
