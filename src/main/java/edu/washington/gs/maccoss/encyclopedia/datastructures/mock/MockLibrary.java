package edu.washington.gs.maccoss.encyclopedia.datastructures.mock;

import java.io.IOException;
import java.nio.file.Path;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;

/**
 * use for testing, does not SQRT!
 * @author searleb
 *
 */
public class MockLibrary implements LibraryInterface {
	private final LibraryEntry[] entries;
		
	public MockLibrary(LibraryEntry[] entries) {
		this.entries=entries;
	}
	
	@Override
	public String getName() {
		return "Testing library of "+entries.length+" entries";
	}
	public ArrayList<LibraryEntry> getAllEntries(boolean sqrt, AminoAcidConstants aaConstants) throws IOException, SQLException, DataFormatException {
		return new ArrayList<LibraryEntry>(Arrays.asList(entries));
	}
	
	public HashMap<String, String> getAccessions(Collection<String> peptideSeqs) throws IOException, SQLException, DataFormatException {
		HashSet<String> set=new HashSet<>(peptideSeqs);
		HashMap<String, String> accessions=new HashMap<>();
		for (LibraryEntry entry : entries) {
			if (set.contains(entry.getPeptideSeq())) {
				accessions.put(entry.getPeptideSeq(), PSMData.accessionsToString(entry.getAccessions()));
			}
		}
		return accessions;
	}
	
	public ArrayList<LibraryEntry> getUnlinkedEntries(Range precursorMz, boolean sqrt, AminoAcidConstants aaConstants) throws IOException, SQLException, DataFormatException {
		return getEntries(precursorMz, sqrt, aaConstants);
	}

	@Override
	public ArrayList<LibraryEntry> getEntries(String peptideModSeq, byte charge, boolean sqrt) throws IOException, SQLException, DataFormatException {
		ArrayList<LibraryEntry> returnables=new ArrayList<LibraryEntry>();
		for (LibraryEntry entry : entries) {
			if (peptideModSeq.equals(entry.getPeptideModSeq())&&charge==entry.getPrecursorCharge()) {
				returnables.add(entry);
			}
		}
		return returnables;
	}

	@Override
	public ArrayList<LibraryEntry> getEntries(Range precursorMz, boolean sqrt, AminoAcidConstants aaConstants) throws IOException, SQLException, DataFormatException {
		ArrayList<LibraryEntry> returnables=new ArrayList<LibraryEntry>();
		for (LibraryEntry entry : entries) {
			if (precursorMz.contains((float)entry.getPrecursorMZ())) {
				returnables.add(entry);
			}
		}
		return returnables;
	}
	
	public Range getMinMaxMZ() throws IOException, SQLException {
		double min=Double.MAX_VALUE;
		double max=0.0;
		for (LibraryEntry entry : entries) {
			if (entry.getPrecursorMZ()>max) max=entry.getPrecursorMZ();
			if (entry.getPrecursorMZ()<min) min=entry.getPrecursorMZ();
		}
		return new Range((float)min, (float)max);
	}
	
	@Override
	public Optional<Path> getSource(SearchParameters parameters) {
		return Optional.empty();
	}
	
	@Override
	public List<Path> getSourceFiles() throws IOException, SQLException {
		return new ArrayList<>();
	}
}
