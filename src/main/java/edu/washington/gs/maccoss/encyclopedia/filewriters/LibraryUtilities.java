package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class LibraryUtilities {
	public static void subsetLibrary(final File saveFile, final float rtMinSec, final float rtMaxSec, final HashSet<String> targets, LibraryFile library)
			throws IOException, SQLException, DataFormatException {
		LibraryFile saveLibrary=new LibraryFile();
		saveLibrary.openFile();
		
		ArrayList<LibraryEntry> toWrite=new ArrayList<>();
		for (LibraryEntry entry : library.getAllEntries(false, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()))) {
			if (rtMinSec<=entry.getRetentionTime()&&rtMaxSec>=entry.getRetentionTime()) {
				if (targets.size()==0||(targets.contains(entry.getPeptideSeq())||targets.contains(entry.getPeptideModSeq()))) {
					toWrite.add(entry);
				}
			}
		}
		Logger.logLine("Found "+toWrite.size()+" peptides from "+targets.size()+" target sequences. Writing to ["+saveFile.getAbsolutePath()+"]...");
		
		saveLibrary.dropIndices();
		saveLibrary.addEntries(toWrite);
		saveLibrary.addProteinsFromEntries(toWrite);
		saveLibrary.createIndices();
		saveLibrary.saveAsFile(saveFile);
		
		saveLibrary.close();
	}
}
