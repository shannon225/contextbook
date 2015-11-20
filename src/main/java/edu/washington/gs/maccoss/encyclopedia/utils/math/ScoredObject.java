package edu.washington.gs.maccoss.encyclopedia.utils.math;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;

public class ScoredObject<T extends Comparable<T>> extends Pair<Float, T> implements Comparable<ScoredObject<T>> {

	public ScoredObject(Float x, T y) {
		super(x, y);
	}

	@Override
	public int compareTo(ScoredObject<T> o) {
		if (o==null) return 1;
		int c=Float.compare(x, o.x);
		if (c!=0) return c;
		
		return y.compareTo(o.y);
	}
}
