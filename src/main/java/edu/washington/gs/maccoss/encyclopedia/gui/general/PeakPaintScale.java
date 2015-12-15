package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.Color;
import java.awt.Paint;
import java.io.Serializable;

import org.jfree.chart.HashUtilities;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.util.PublicCloneable;

public class PeakPaintScale implements PaintScale, PublicCloneable, Serializable {
	private static final long serialVersionUID=1L;
	
	final private double lowerBound;
	final private double upperBound;
	final private double range;
	
	final Color[] colors=new Color[256];

	public PeakPaintScale(double lowerBound, double upperBound) {
		this(lowerBound, upperBound, 255);
	}

	public PeakPaintScale(double lowerBound, double upperBound, int alpha) {
		if (lowerBound>=upperBound) {
			throw new IllegalArgumentException("Requires lowerBound < upperBound.");
		}
		if (alpha<0||alpha>255) {
			throw new IllegalArgumentException("Requires alpha in the range 0 to 255.");
		}
		this.lowerBound=lowerBound;
		this.upperBound=upperBound;
		this.range=upperBound-lowerBound;
		
		for (int i=0; i<colors.length; i++) {
			colors[255-i]=new Color(i, i, i, alpha);
		}
	}

	@Override
	public double getLowerBound() {
		return this.lowerBound;
	}

	@Override
	public double getUpperBound() {
		return this.upperBound;
	}

	@Override
	public Paint getPaint(double value) {
		double v=Math.min(Math.max(value, lowerBound), upperBound);
		int g=(int)((v-lowerBound)/(range)*255.0);
		
		return colors[g];
	}

	@Override
	public boolean equals(Object obj) {
		if (obj==this) {
			return true;
		}
		if (!(obj instanceof PeakPaintScale)) {
			return false;
		}
		PeakPaintScale that=(PeakPaintScale)obj;
		if (this.lowerBound!=that.lowerBound) {
			return false;
		}
		if (this.upperBound!=that.upperBound) {
			return false;
		}
		return true;
	}

	@Override
	public int hashCode() {
		int hash=7;
		hash=HashUtilities.hashCode(hash, this.lowerBound);
		hash=HashUtilities.hashCode(hash, this.upperBound);
		return hash;
	}

	@Override
	public Object clone() throws CloneNotSupportedException {
		return super.clone();
	}

}
