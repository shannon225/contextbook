package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public class PeptideModification {

	public static final String NO_MODIFICATION_NAME = "none";
	public static PeptideModification phosphorylation = new PeptideModification("Phosphorylation (STY)", "Phosphorylation", 79.966331, new double[]{97.976896, 97.976896, 0.0}, new char[]{'S', 'T', 'Y'});
	static PeptideModification polymorphism           = new PeptideModification("Nucleotide Polymorphism", "Polymorphism", 0.0, new double["ACDEFGHIKLMNPQRSTVWY".length()], "ACDEFGHIKLMNPQRSTVWY".toCharArray());

	private final String name;
	private final String shortname;
	private final double mass;
	private final double[] neutralLoss; // one mass per modifiableAA
	private final int nominalMass;
	private final char[] modifiableAAs;

	public PeptideModification(String name, String shortname, double mass, double[] neutralLoss, char[] modifiableAAs) {
		this.name=name;
		this.shortname=shortname;
		this.mass=mass;
		this.neutralLoss=neutralLoss;
		this.nominalMass=(int)Math.round(mass);
		this.modifiableAAs=modifiableAAs;
	}
	@Override
	public String toString() {
		return name+" ["+(mass>0?"+":"")+mass+"]";
	}
	public String toMassString() {
		return "["+(mass>0?"+":"")+mass+"]";
	}
	public double getMass() {
		return mass;
	}
	public char[] getModifiableAAs() {
		return modifiableAAs;
	}
	public String getName() {
		return name;
	}
	public String getShortname() {
		return shortname;
	}
	public int getNominalMass() {
		return nominalMass;
	}
	
	private static final MassTolerance tolerance=new MassTolerance(1.0); // 1 ppm is about the accuracy of floats 

	public boolean isModifiable(char aa) {
		for (int i=0; i<modifiableAAs.length; i++) {
			if (modifiableAAs[i]==aa) {
				return true;
			}
		}
		return false;
	}
	public boolean isModificationMass(char aa, double modificationMass) {
		boolean ok=false;
		for (int i=0; i<modifiableAAs.length; i++) {
			if (modifiableAAs[i]==aa) {
				ok=true;
			}
		}
		if (ok) {
			if (tolerance.equals(mass, modificationMass)) {
				return true;
			} else if (tolerance.equals((double) nominalMass, modificationMass)) {
				return true;
			}
		}
		return false;
		
	}

	public double getNeutralLoss(char aa) {
		for (int i=0; i<modifiableAAs.length; i++) {
			if (modifiableAAs[i]==aa) {
				if (i>=neutralLoss.length) {
					return 0.0;
				} else {
					return neutralLoss[i];
				}
			}
		}
		return 0.0;
	}

}
