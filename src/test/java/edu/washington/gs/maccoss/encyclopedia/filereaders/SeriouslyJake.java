package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;

public class SeriouslyJake {
	public static void main(String[] args) throws Exception {
		checkYourData();
	}
	public static void checkYourData() throws Exception {
		LibraryFile library=new LibraryFile();
		library.openFile(new File("/Users/searleb/Documents/projects/p100/PC3_DDA_UnambigPhosphoLibrary.nonredundant.AllTables.dlib"));
		
		ArrayList<LibraryEntry> entries=library.getAllEntries(false, new AminoAcidConstants());
		HashMap<String, LibraryEntry> entryMap=new HashMap<>();
		for (LibraryEntry entry : entries) {
			String key=entry.getPeptideModSeq()+"+"+entry.getPrecursorCharge()+"+"+entry.getTIC();
			entryMap.put(key, entry);
		}
		
		LibraryFile newLibrary=new LibraryFile();
		newLibrary.openFile();
		newLibrary.addEntries(new ArrayList<LibraryEntry>(entryMap.values()));
		newLibrary.saveAsFile(new File("/Users/searleb/Documents/projects/p100/cleaned_PC3_DDA_UnambigPhosphoLibrary.nonredundant.AllTables.dlib"));
		System.out.println("Wrote "+entryMap.size()+" unique entries of "+entries.size()+" total entries");
	}
}
