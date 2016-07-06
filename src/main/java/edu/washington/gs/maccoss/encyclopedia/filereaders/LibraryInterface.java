package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;

public interface LibraryInterface {
	ArrayList<LibraryEntry> getEntries(String peptideModSeq, byte charge, boolean sqrt) throws IOException, SQLException, DataFormatException;
	ArrayList<LibraryEntry> getEntries(Range precursorMz, boolean sqrt) throws IOException, SQLException, DataFormatException;
	ArrayList<LibraryEntry> getAllEntries(boolean sqrt) throws IOException, SQLException, DataFormatException;
	Range getMinMaxMZ() throws IOException, SQLException;
	String getName();
	Optional<StripeFileInterface> getSource(SearchParameters parameters);
}