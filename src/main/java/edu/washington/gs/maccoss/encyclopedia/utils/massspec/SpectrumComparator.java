package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.Comparator;

public class SpectrumComparator implements Comparator<Spectrum> {
	public static final byte compareWithName=0;
	public static final byte compareWithMZ=1;
	public static final byte compareWithRT=2;
	
	private final byte compareMethod;

	public SpectrumComparator(byte compareMethod) {
		this.compareMethod=compareMethod;
	}
	
	@Override
	public int compare(Spectrum o1, Spectrum o2) {
		if (o1==null&&o2==null) return 0;
		if (o1==null) return -1;
		if (o2==null) return 1;
		
		if (compareMethod==compareWithName) {
			return o1.getSpectrumName().compareTo(o2.getSpectrumName()); 
		} else if (compareMethod==compareWithMZ) {
			return Double.compare(o1.getPrecursorMZ(), o2.getPrecursorMZ());
		} else {
			return Float.compare(o1.getScanStartTime(), o2.getScanStartTime());
		}
	}
}
