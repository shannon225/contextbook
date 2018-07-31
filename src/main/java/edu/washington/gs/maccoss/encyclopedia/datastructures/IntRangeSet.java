package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;

public class IntRangeSet {
	private final ArrayList<IntRange> ranges=new ArrayList<>();

	public IntRangeSet() {
	}

	public void addRange(IntRange range) {
		ranges.add(range);
	}
	
	public boolean contains(int index) {
		for (IntRange range : ranges) {
			if (range.contains(index)) return true;
		}
		return false;
	}
}
