package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.awt.Color;

import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;

public class PrecursorIon implements Comparable<PrecursorIon>, Ion {
	private static final MassTolerance tolerance = new MassTolerance(0.1); // 1 ppm is about the accuracy of floats

	private final String annotation;
	private final double mass;
	private final byte charge;
	private final byte isotope;

	public PrecursorIon(String annotation, double mass, byte charge) {
		this.annotation = annotation;
		this.mass = mass;
		this.charge=charge;
		this.isotope=0;
	}

	public PrecursorIon(String annotation, double mass, byte charge, byte isotope) {
		this.annotation = annotation;
		this.mass = mass;
		this.charge=charge;
		this.isotope=isotope;
	}

	public double getMass() {
		return mass;
	}

	@Override
	public String getName() {
		return annotation;
	}
	
	public byte getCharge() {
		return charge;
	}
	
	public byte getIsotope() {
		return isotope;
	}
	
	@Override
	public byte getIndex() {
		return isotope;
	}
	
	@Override
	public IonType getType() {
		return IonType.precursor;
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
