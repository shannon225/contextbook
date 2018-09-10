package edu.washington.gs.maccoss.encyclopedia.utils;

public class ComparablePair<X extends Comparable<? super X>, Y extends Comparable<? super Y>> extends Pair<X, Y>
		implements Comparable<ComparablePair<X, Y>> {

	public ComparablePair(X x, Y y) {
		super(x, y);
	}

	@Override
	public int compareTo(ComparablePair<X, Y> o) {
		int firstCompareTo = x.compareTo(o.x);

		if (firstCompareTo == 0) {
			return y.compareTo(o.y);
		}

		return firstCompareTo;
	}
	
	@Override
	public int hashCode() {
		return x.hashCode()+16807*y.hashCode();
	}
	
	@SuppressWarnings("rawtypes")
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof ComparablePair) {
			if (!x.equals(((ComparablePair)obj).x)) {
				return false;
			}
			if (!y.equals(((ComparablePair)obj).y)) {
				return false;
			}
			return true;
		}
		return false;
	}
}