package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class LibraryUtilitiesTest {
	public static void main2(String[] args) throws Exception {

		File twoDLC=new File("/Users/searleb/Downloads/msms.dlib");

		LibraryFile library=new LibraryFile();
		library.openFile(twoDLC);

		float[] timeBoundaries=new float[] {0f, 2f, 32f, 62f, 92f, 122f, 152f, 182f, 212f, 242f, 272f, 302f, 332f, 362f, 392f, 1000f};
		for (int i=1; i<timeBoundaries.length; i++) {
			File saveFile=new File(twoDLC.getParentFile(), "msms_"+Math.round(timeBoundaries[i-1])+"to"+Math.round(timeBoundaries[i])+".dlib");
			System.out.println("writing "+saveFile.getName());
			LibraryUtilities.subsetLibrary(saveFile, timeBoundaries[i-1]*60f, timeBoundaries[i]*60f, new HashSet<>(), library);
		}

		library.close();
	}
	
	public static void main(String[] args) throws Exception {
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
}
