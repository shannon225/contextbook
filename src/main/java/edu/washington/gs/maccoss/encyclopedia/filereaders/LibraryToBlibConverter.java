package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import gnu.trove.map.hash.TCharDoubleHashMap;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.zip.DataFormatException;

public class LibraryToBlibConverter {

	public static void convert(File elibFile, File blibFile) throws EncyclopediaException {
		try {
			LibraryFile library = new LibraryFile();
			library.openFile(elibFile);
			final AminoAcidConstants aaConstants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
			final ArrayList<LibraryEntry> allEntries = library.getAllEntries(false, aaConstants);
			Logger.logLine("Found " + allEntries.size() + " entries from " + elibFile.getName() + ". Writing to [" + blibFile.getAbsolutePath() + "]...");

			BlibFile blib = new BlibFile();
			blib.openFile();
			blib.setUserFile(blibFile);
			blib.dropIndices();
			int[] counterTotals = new int[]{0, 0, 0};

			blib.addLibrary(allEntries, library.getName(), aaConstants, "ELIB conversion", counterTotals[0], counterTotals[1], counterTotals[2]);
			blib.createIndices();
			blib.saveFile();
			blib.close();
			library.close();
			Logger.logLine("Finished reading " + blibFile.getName());
		} catch (IOException | SQLException | DataFormatException e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
			throw new EncyclopediaException("Encountered fatal error converting library", e);
		}
	}
}
