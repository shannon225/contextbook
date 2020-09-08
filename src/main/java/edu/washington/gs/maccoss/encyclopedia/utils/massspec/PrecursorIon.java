package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.awt.Color;

import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;

public class PrecursorIon implements Comparable<PrecursorIon>, Ion {
	private static final MassTolerance tolerance = new MassTolerance(0.1); // 1 ppm is about the accuracy of floats

	private final String annotation;
	private final double mass;

	public PrecursorIon(String annotation, double mass) {
		this.annotation = annotation;
		this.mass = mass;
	}

	public double getMass() {
		return mass;
	}

	@Override
	public String getName() {
		return annotation;
	}

	@Override
	public Color getColor() {
		return RandomGenerator.randomColor(getName().hashCode());
	}

	@Override
	public int hashCode() {
		// Note that equal objects will always have identical masses (see below)
		return (int)(mass*100.0);
	}

	@Override
	public int compareTo(PrecursorIon o) {
		if (o == null) {
			return 1;
		}
		int c = tolerance.compareTo(mass, o.mass);
		if (c != 0)
			return c;
		return annotation.compareTo(o.annotation);
	}
	
	@Override
	public String toString() {
		return getName();
	}

}
