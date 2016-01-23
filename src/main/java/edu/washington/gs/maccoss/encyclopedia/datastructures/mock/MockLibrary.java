package edu.washington.gs.maccoss.encyclopedia.datastructures.mock;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import jdk.nashorn.internal.ir.annotations.Immutable;

/**
 * use for testing, does not SQRT!
 * @author searleb
 *
 */
public class MockLibrary implements LibraryInterface {
	@Immutable
	private final LibraryEntry[] entries;
		
	public MockLibrary(LibraryEntry[] entries) {
		this.entries=entries;
	}

	@Override
	public ArrayList<LibraryEntry> getEntries(String peptideModSeq, boolean sqrt) throws IOException, SQLException, DataFormatException {
		ArrayList<LibraryEntry> returnables=new ArrayList<LibraryEntry>();
		for (LibraryEntry entry : entries) {
			if (peptideModSeq.equals(entry.getPeptideModSeq())) {
				returnables.add(entry);
			}
		}
		return returnables;
	}

	@Override
	public ArrayList<LibraryEntry> getEntries(Range precursorMz, boolean sqrt) throws IOException, SQLException, DataFormatException {
		ArrayList<LibraryEntry> returnables=new ArrayList<LibraryEntry>();
		for (LibraryEntry entry : entries) {
			if (precursorMz.contains((float)entry.getPrecursorMZ())) {
				returnables.add(entry);
			}
		}
		return returnables;
	}

}
