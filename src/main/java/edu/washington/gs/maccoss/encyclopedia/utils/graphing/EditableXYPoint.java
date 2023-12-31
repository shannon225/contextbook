package edu.washington.gs.maccoss.encyclopedia.utils.graphing;

/**
 * Unstable (but safe-ish) if used in a non-single threaded way!
 * If values are null then getX and getY return NaN
 * @author searleb
 *
 */
public class EditableXYPoint implements PointInterface {
	private Double x;
	private Double y;
	
	public EditableXYPoint() {
		x=null;
		y=null;
	}

	public EditableXYPoint(double x, double y) {
		this.x=x;
		this.y=y;
	}
	
	@Override
	public String toString() {
		return x+","+y;
	}

	@Override
	public double getX() {
		if (x==null) return Double.NaN;
		return x;
	}

	@Override
	public double getY() {
		if (y==null) return Double.NaN;
		return y;
	}
	
	public void setX(Double x) {
		this.x=x;
	}

	public void setY(Double y) {
		this.y=y;
	}

	/**
	 * compares on X first then on Y
	 */
	@Override
	public int compareTo(PointInterface o) {
		if (o==null) return 1;
		int c=Double.compare(getX(), o.getX());
		if (c!=0) return c;
		c=Double.compare(getY(), o.getY());
		return c;
	}
	
	@Override
	public int hashCode() {
		return Double.hashCode(getX())+Double.hashCode(getY());
	}
	
	@Override
	public boolean equals(Object obj) {
		return compareTo((PointInterface)obj)==0;
	}
}
