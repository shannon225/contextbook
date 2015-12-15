package edu.washington.gs.maccoss.encyclopedia.utils.math;

public class ScoredIndex implements Comparable<ScoredIndex> {
	public final float x;
	public final int y;

	public ScoredIndex(float x, int y) {
		this.x=x;
		this.y=y;
	}

	@Override
	public int compareTo(ScoredIndex o) {
		if (o==null) return 1;
		int c=Float.compare(x, o.x);
		if (c!=0) return c;

		return Integer.compare(y, o.y);
	}
}
