package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;

public interface LibraryInterface {
	ArrayList<LibraryEntry> getEntries(String peptideModSeq, boolean sqrt) throws IOException, SQLException, DataFormatException;
	ArrayList<LibraryEntry> getEntries(Range precursorMz, boolean sqrt) throws IOException, SQLException, DataFormatException;

}