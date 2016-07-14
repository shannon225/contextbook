package edu.washington.gs.maccoss.encyclopedia.datastructures.mock;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
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
	public String getName() {
		return "Testing library of "+entries.length+" entries";
	}
	public ArrayList<LibraryEntry> getAllEntries(boolean sqrt) throws IOException, SQLException, DataFormatException {
		return new ArrayList<LibraryEntry>(Arrays.asList(entries));
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
	public HashMap<PeptidePrecursor, ArrayList<LibraryEntry>> getEntries(ArrayList<PeptidePrecursor> entries, boolean sqrt) throws IOException, SQLException, DataFormatException {
		HashMap<PeptidePrecursor, ArrayList<LibraryEntry>> map=new HashMap<PeptidePrecursor, ArrayList<LibraryEntry>>();
		for (PeptidePrecursor peptidePrecursor : entries) {
			ArrayList<LibraryEntry> results=getEntries(peptidePrecursor.getPeptideModSeq(),  peptidePrecursor.getPrecursorCharge(), sqrt);
			map.put(peptidePrecursor, results);
		}
		return map;
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
	
	@Override
	public Range getMinMaxMZ() throws IOException, SQLException {
		double min=Double.MAX_VALUE;
		double max=0.0;
		for (LibraryEntry entry : entries) {
			if (entry.getPrecursorMZ()>max) max=entry.getPrecursorMZ();
			if (entry.getPrecursorMZ()<min) min=entry.getPrecursorMZ();
		}
		return new Range((float)min, (float)max);
	}
	
	public Optional<StripeFileInterface> getSource(SearchParameters parameters) {
		return Optional.empty();
	}
}
