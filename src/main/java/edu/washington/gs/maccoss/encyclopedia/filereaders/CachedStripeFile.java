package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;

public class CachedStripeFile implements StripeFileInterface {
	private final File userFile;
	private final HashMap<Range, Float> ranges;
	private final ArrayList<PrecursorScan> precursors;
	private final HashMap<Range, ArrayList<Stripe>> stripes;
	
	public CachedStripeFile(File userFile, HashMap<Range, Float> ranges, ArrayList<PrecursorScan> precursors, HashMap<Range, ArrayList<Stripe>> stripes) {
		this.userFile=userFile;
		this.ranges=ranges;
		this.precursors=precursors;
		this.stripes=stripes;
	}

	@Override
	public HashMap<Range, Float> getRanges() {
		return ranges;
	}

	@Override
	public void openFile(File userFile) throws IOException, SQLException {
	}

	@Override
	public void openFile() throws IOException, SQLException {
	}

	@Override
	public ArrayList<PrecursorScan> getPrecursors(float minRT, float maxRT) throws IOException, SQLException, DataFormatException {
		ArrayList<PrecursorScan> subset=new ArrayList<PrecursorScan>();
		for (PrecursorScan scan : precursors) {
			if (scan.getScanStartTime()>=minRT&&scan.getScanStartTime()<=maxRT) {
				subset.add(scan);
			}
		}
		return subset;
	}

	@Override
	public ArrayList<Stripe> getStripes(double targetMz, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException {
		for (Entry<Range, ArrayList<Stripe>> entry : stripes.entrySet()) {
			if (entry.getKey().contains((float)targetMz)) {
				ArrayList<Stripe> subset=new ArrayList<Stripe>();
				for (Stripe scan : entry.getValue()) {
					if (scan.getScanStartTime()>=minRT&&scan.getScanStartTime()<=maxRT) {
						if (sqrt) {
							subset.add(scan.sqrt());
						} else {
							subset.add(scan);
						}
					}
				}
				return subset;
			}
		}
		
		return new ArrayList<Stripe>();
	}
	
	@Override
	public ArrayList<Stripe> getStripes(Range targetMzRange, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException {
		for (Entry<Range, ArrayList<Stripe>> entry : stripes.entrySet()) {
			if (targetMzRange.contains(entry.getKey().getMiddle())) {
				ArrayList<Stripe> subset=new ArrayList<Stripe>();
				for (Stripe scan : entry.getValue()) {
					if (scan.getScanStartTime()>=minRT&&scan.getScanStartTime()<=maxRT) {
						if (sqrt) {
							subset.add(scan.sqrt());
						} else {
							subset.add(scan);
						}
					}
				}
				return subset;
			}
		}
		
		return new ArrayList<Stripe>();
	}

	@Override
	public void close() {
		ranges.clear();
		precursors.clear();
		stripes.clear();
	}

	@Override
	public File getFile() {
		return userFile;
	}

}
