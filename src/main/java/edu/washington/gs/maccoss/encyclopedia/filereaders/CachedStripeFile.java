package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;

public class CachedStripeFile implements StripeFileInterface {
	private final File userFile;
	private final Map<String, String> metadata;
	private final Map<Range, WindowData> ranges;
	private final List<PrecursorScan> precursors;
	private final Map<Range, ? extends List<FragmentScan>> stripes;
	private final float tic;
	private final float gradientLength;
	
	public CachedStripeFile(File userFile, Map<Range, WindowData> ranges,  Map<String, String> metadata, List<PrecursorScan> precursors, Map<Range, ? extends List<FragmentScan>> stripes) {
		this.userFile=userFile;
		this.metadata=metadata;
		this.ranges=ranges;
		this.precursors=precursors;
		this.stripes=stripes;
		
		float sum=0.0f;
		for (PrecursorScan precursorScan : precursors) {
			sum+=precursorScan.getTIC();
		}
		tic=sum;
		
		float maxRT=0.0f;
		for (List<FragmentScan> stripe : stripes.values()) {
			for (FragmentScan scan : stripe) {
				if (scan.getScanStartTime()>maxRT) {
					maxRT=scan.getScanStartTime();
				}
			}
		}
		gradientLength=maxRT;
	}
	
	@Override
	public Map<String, String> getMetadata() throws IOException, SQLException {
		return metadata;
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
	public Map<Range, WindowData> getRanges() {
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
	public ArrayList<FragmentScan> getStripes(double targetMz, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException {
		for (Entry<Range, ? extends List<FragmentScan>> entry : stripes.entrySet()) {
			if (entry.getKey().contains((float)targetMz)) {
				ArrayList<FragmentScan> subset=new ArrayList<FragmentScan>();
				for (FragmentScan scan : entry.getValue()) {
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
		
		return new ArrayList<FragmentScan>();
	}
	
	@Override
	public ArrayList<FragmentScan> getStripes(Range targetMzRange, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException {
		for (Entry<Range, ? extends List<FragmentScan>> entry : stripes.entrySet()) {
			if (targetMzRange.contains(entry.getKey().getMiddle())) {
				ArrayList<FragmentScan> subset=new ArrayList<FragmentScan>();
				for (FragmentScan scan : entry.getValue()) {
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
		
		return new ArrayList<FragmentScan>();
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
