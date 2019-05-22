package edu.washington.gs.maccoss.encyclopedia.utils.graphing;


public class XYZPoint implements PointInterface {
	public final double x;
	public final double y;
	public final double z;

	public XYZPoint(double x, double y, double z) {
		this.x=x;
		this.y=y;
		this.z=z;
	}

	@Override
	public double getX() {
		return x;
	}

	@Override
	public double getY() {
		return y;
	}
	
	public double getZ() {
		return z;
	}

	/**
	 * compares on X first then on Y. Only try Z if other PointInterface is an XYZPoint
	 */
	@Override
	public int compareTo(PointInterface o) {
		if (o==null) return 1;
		if (x>o.getX()) return 1;
		if (x<o.getX()) return -1;
		if (y>o.getY()) return 1;
		if (y<o.getY()) return -1;
		if (o instanceof XYZPoint) {
			if (z>((XYZPoint)o).getZ()) return 1;
			if (z<((XYZPoint)o).getZ()) return -1;
		}
		return 0;
	}
}
