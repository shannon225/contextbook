package edu.washington.gs.maccoss.encyclopedia.utils.math;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;

public class IndexedObject<T extends Comparable<T>> extends Pair<Integer, T> implements Comparable<IndexedObject<T>> {

	public IndexedObject(Integer x, T y) {
		super(x, y);
	}

	@Override
	public int compareTo(IndexedObject<T> o) {
		if (o==null) return 1;
		int c=Integer.compare(x, o.x);
		if (c!=0) return c;
		
		return y.compareTo(o.y);
	}
}
