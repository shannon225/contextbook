package edu.washington.gs.maccoss.encyclopedia.utils.graphing;

import jdk.nashorn.internal.ir.annotations.Immutable;

@Immutable
public class XYPoint implements PointInterface {
	public final double x;
	public final double y;

	public XYPoint(double x, double y) {
		this.x=x;
		this.y=y;
	}

	@Override
	public double getX() {
		return x;
	}

	@Override
	public double getY() {
		return y;
	}

	/**
	 * compares on X first then on Y
	 */
	@Override
	public int compareTo(PointInterface o) {
		if (o==null) return 1;
		if (x>o.getX()) return 1;
		if (x<o.getX()) return -1;
		if (y>o.getY()) return 1;
		if (y<o.getY()) return -1;
		return 0;
	}
}
