package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.IOException;
import java.nio.file.Path;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.zip.DataFormatException;

import com.google.common.util.concurrent.AtomicDouble;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;

public class OrderedCachedLibraryInterface implements LibraryInterface {
	private final LibraryInterface source;
	
	// FIXME TODO!
	// tries to cache entries keeping track of our current cache, assuming requests are always increasing
	private volatile AtomicDouble currentMinimum=new AtomicDouble(0.0);
	private volatile AtomicDouble currentMaximum=new AtomicDouble(0.0);
	
	public OrderedCachedLibraryInterface(LibraryInterface source) {
		this.source=source;
	}
	@Override
	public String getName() {
		return source.getName();
	}
	@Override
	public Optional<Path> getSource(SearchParameters parameters) {
		return source.getSource(parameters);
	}
	@Override
	public List<Path> getSourceFiles() throws IOException, SQLException {
		return source.getSourceFiles();
	}

	
	@Override
	public ArrayList<LibraryEntry> getEntries(String peptideModSeq, byte charge, boolean sqrt)
			throws IOException, SQLException, DataFormatException {
		// TODO Auto-generated method stub
		return source.getEntries(peptideModSeq, charge, sqrt);
	}
	@Override
	public ArrayList<LibraryEntry> getEntries(Range precursorMz, boolean sqrt, AminoAcidConstants aaConstants)
			throws IOException, SQLException, DataFormatException {
		// TODO Auto-generated method stub
		return source.getEntries(precursorMz, sqrt, aaConstants);
	}
	@Override
	public ArrayList<LibraryEntry> getAllEntries(boolean sqrt, AminoAcidConstants aaConstants)
			throws IOException, SQLException, DataFormatException {
		// TODO Auto-generated method stub
		return source.getAllEntries(sqrt, aaConstants);
	}
	
	
	
}
