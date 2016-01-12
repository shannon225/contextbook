package edu.washington.gs.maccoss.encyclopedia.utils;

public class Pair<X, Y> {
	public final X x;
	public final Y y;

	public Pair(X x, Y y) {
		this.x = x;
		this.y = y;
	}
	
	@Override
	public String toString() {
		return x+" and "+y;
	}
}