package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.SearchTestSupport;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import junit.framework.TestCase;

public class LibraryEntryCleanerTest extends TestCase {
	public void testFilterIons() throws Exception {
		LibraryFile library = (LibraryFile)BlibToLibraryConverter.getFile(SearchTestSupport.getNewTestLibraryFile());
		ArrayList<LibraryEntry> newList=LibraryEntryCleaner.filterIons(library.getAllEntries(false, new AminoAcidConstants()), 5);
		for (LibraryEntry entry : newList) {
			assertTrue(entry.getMassArray().length<=5);
		}
	}
}
