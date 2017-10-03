package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class CachedStripeFile implements StripeFileInterface {
	private final File userFile;
	private final Map<Range, Float> ranges;
	private final List<PrecursorScan> precursors;
	private final Map<Range, ? extends List<Stripe>> stripes;
	private final float tic;
	private final float gradientLength;
	
	public CachedStripeFile(File userFile, Map<Range, Float> ranges, List<PrecursorScan> precursors, Map<Range, ? extends List<Stripe>> stripes) {
		this.userFile=userFile;
		this.ranges=ranges;
		this.precursors=precursors;
		this.stripes=stripes;
		
		float sum=0.0f;
		for (PrecursorScan precursorScan : precursors) {
			sum+=General.sum(precursorScan.getIntensityArray());
		}
		tic=sum;
		
		float maxRT=0.0f;
		for (List<Stripe> stripe : stripes.values()) {
			for (Stripe scan : stripe) {
				if (scan.getScanStartTime()>maxRT) {
					maxRT=scan.getScanStartTime();
				}
			}
		}
		gradientLength=maxRT;
	}
	
	@Override
	public float getTIC() throws IOException, SQLException {
		return tic;
	}
	
	@Override
	public float getGradientLength() throws IOException, SQLException {
		return gradientLength;
	}

	@Override
	public Map<Range, Float> getRanges() {
		return ranges;
	}

	@Override
	public void openFile(File userFile) throws IOException, SQLException {
	}
	
	@Override
	public String getOriginalFileName() {
		return userFile.getName();
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
		for (Entry<Range, ? extends List<Stripe>> entry : stripes.entrySet()) {
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
		for (Entry<Range, ? extends List<Stripe>> entry : stripes.entrySet()) {
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
	public boolean isOpen() {
		return ranges.size()>0;
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
