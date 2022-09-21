package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.IOException;
import java.nio.file.Path;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class OrderedCachedLibraryInterface implements LibraryInterface {
	private final LibraryInterface source;
	
	// tries to cache entries keeping track of our current cache, assuming requests are always increasing
	private volatile double currentMaximum=-Double.MIN_VALUE;
	private volatile boolean gotAllValues=false;
	private final ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
	
	private final Comparator<LibraryEntry> mzComparator=new Comparator<LibraryEntry>() {
		@Override
		public int compare(LibraryEntry o1, LibraryEntry o2) {
			if (o1==null&&o2==null) return 0;
			if (o1==null) return -1;
			if (o2==null) return 1;
			return Double.compare(o1.getPrecursorMZ(), o2.getPrecursorMZ());
		}
	};
	
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
	public synchronized ArrayList<LibraryEntry> getEntries(String peptideModSeq, byte charge, boolean sqrt)
			throws IOException, SQLException, DataFormatException {
		// don't bother caching single entries
		return source.getEntries(peptideModSeq, charge, sqrt);
	}
	
	@Override
	public synchronized ArrayList<LibraryEntry> getAllEntries(boolean sqrt, AminoAcidConstants aaConstants)
			throws IOException, SQLException, DataFormatException {
		if (gotAllValues) return entries;
		entries.clear();
		entries.addAll(source.getAllEntries(sqrt, aaConstants));
		Collections.sort(entries, mzComparator);
		gotAllValues=true;
		return entries;
	}
	
	@Override
	public synchronized ArrayList<LibraryEntry> getEntries(Range precursorMz, boolean sqrt, AminoAcidConstants aaConstants)
			throws IOException, SQLException, DataFormatException {
		Logger.logLine("Ordered entry cache currently tracking up to "+(Math.round(currentMaximum*10.0f)/10.0f)+", adjusting to "+precursorMz.toString());
		
		gotAllValues=false;
		ArrayList<LibraryEntry> newList=new ArrayList<LibraryEntry>();
		for (LibraryEntry entry : entries) {
			if (precursorMz.contains(entry.getPrecursorMZ())) newList.add(entry);
		}
		entries.clear();
		entries.addAll(newList);
		
		double nextMax=precursorMz.getStop();
		if (currentMaximum<nextMax) {
			Range shorterRange=new Range(currentMaximum, nextMax);
			ArrayList<LibraryEntry> nextEntries=source.getEntries(shorterRange, sqrt, aaConstants);
			Collections.sort(nextEntries, mzComparator);
			entries.addAll(nextEntries);
		}
		currentMaximum=precursorMz.getStop();
		
		return entries;
	}
	
	
	
}
