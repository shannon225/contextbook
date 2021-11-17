package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.IOException;
import java.nio.file.Path;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;

public interface LibraryInterface {
	ArrayList<LibraryEntry> getEntries(String peptideModSeq, byte charge, boolean sqrt) throws IOException, SQLException, DataFormatException;
//	HashMap<PeptidePrecursor, ArrayList<LibraryEntry>> getEntries(ArrayList<PeptidePrecursor> entries, boolean sqrt) throws IOException, SQLException, DataFormatException;
	ArrayList<LibraryEntry> getEntries(Range precursorMz, boolean sqrt, AminoAcidConstants aaConstants) throws IOException, SQLException, DataFormatException;
	ArrayList<LibraryEntry> getUnlinkedEntries(Range precursorMz, boolean sqrt, AminoAcidConstants aaConstants) throws IOException, SQLException, DataFormatException;
	ArrayList<LibraryEntry> getAllEntries(boolean sqrt, AminoAcidConstants aaConstants) throws IOException, SQLException, DataFormatException;
	HashMap<String, String> getAccessions(Collection<String> peptideSeqs) throws IOException, SQLException, DataFormatException;
	Range getMinMaxMZ() throws IOException, SQLException;
	String getName();
	Optional<Path> getSource(SearchParameters parameters);
	public List<Path> getSourceFiles() throws IOException, SQLException;
}