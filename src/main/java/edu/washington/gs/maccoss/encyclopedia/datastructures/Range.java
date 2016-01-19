package edu.washington.gs.maccoss.encyclopedia.datastructures;


//@Immutable
public class Range implements Comparable<Range> {
	private final float start, stop;

	public Range(float start, float stop) {
		// ensure that start comes before stop
		if (start<=stop) {
			this.start = start;
			this.stop = stop;
		} else {
			this.start = stop;
			this.stop = start;
		}
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
	
	public float getStart() {
		return start;
	}
	
	public float getStop() {
		return stop;
	}
	
	public float getMiddle() {
		return (start+stop)/2.0f;
	}
	
	public float getRange() {
		return stop-start;
	}
	
	public boolean contains(float value) {
		if (value>=start&&value<=stop) {
			return true;
		}
		return false;
	}
	
	/**
	 * sorts on start location, then on stop location
	 */
	public int compareTo(Range o) {
		if (o==null) return 1;
		if (start>o.start) return 1;
		if (start<o.start) return -1;
		if (stop>o.stop) return 1;
		if (stop<o.stop) return -1;
		return 0;
	}
}
