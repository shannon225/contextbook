package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;

//@Immutable
public class IntRange implements Comparable<Range> {
	private final int start, stop;

	public IntRange(int start, int stop) {
		// ensure that start comes before stop
		if (start<=stop) {
			this.start = start;
			this.stop = stop;
		} else {
			this.start = stop;
			this.stop = start;
		}
	}
	
	public int getLength() {
		return stop-start+1;
	}
	
	@Override
	public String toString() {
		return start+" to "+stop;
	}
	
	@Override
	public int hashCode() {
		return Float.floatToIntBits(start)+16807*Float.floatToIntBits(stop);
	}
	
	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof Range)) return false;
		return compareTo((Range)obj)==0;
	}
	
	public int getStart() {
		return start;
	}
	
	public int getStop() {
		return stop;
	}
	
	public int getMiddle() {
		return (start+stop)/2;
	}
	
	public int getRange() {
		return stop-start;
	}
	
	public boolean contains(float value) {
		if (value>=start&&value<=stop) {
			return true;
		}
		return false;
	}
	
	public ArrayList<IntRange> chunkIntoBins(int binCount) {
		int delta=getRange()/binCount;
		
		ArrayList<IntRange> ranges=new ArrayList<IntRange>();
		int currentMin=start;
		for (int i=0; i<binCount; i++) {
			int currentMax=currentMin+delta;
			ranges.add(new IntRange(currentMin, currentMax));
			currentMin=currentMax;
		}
		return ranges;
	}
	
	/**
	 * sorts on start location, then on stop location
	 */
	public int compareTo(Range o) {
		if (o==null) return 1;
		if (start>o.getStart()) return 1;
		if (start<o.getStart()) return -1;
		if (stop>o.getStop()) return 1;
		if (stop<o.getStop()) return -1;
		return 0;
	}
	
	public float linearInterp(float X, float minY, float maxY) {
		float deltaX=getRange();
		if (deltaX==0) {
			float half=(maxY+minY)/2f;
			if (half<minY) return minY;
			if (half>maxY) return maxY;
			return half;
		}
		float deltaY=maxY-minY;
		if (deltaY==0) {
			return maxY;
		}
		float interp=((deltaY/deltaX)*(X-getStart())+minY);
		if (interp<minY) return minY;
		if (interp>maxY) return maxY;
		return interp;
	}
	
	public int linearInterp(float X, int minY, int maxY) {
		return Math.round(linearInterp(X, (float)minY, (float)maxY));
	}
	
	public float mapBackToRange(float Y, float minY, float maxY) {
		float deltaX=getRange();
		if (deltaX==0) {
			return getStop();
		}
		
		float deltaY=maxY-minY;
		if (deltaY==0) {
			float half=(getStart()+getStop())/2f;
			if (half<getStart()) return getStart();
			if (half>getStop()) return getStop();
			return half;
		}
		float interp=((deltaX/deltaY)*(Y-minY)+getStart());
		if (interp<getStart()) return getStart();
		if (interp>getStop()) return getStop();
		return interp;
	}
}
