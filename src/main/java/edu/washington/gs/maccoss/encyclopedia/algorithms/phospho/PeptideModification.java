package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public class PeptideModification {
	public static final String NO_MODIFICATION_NAME="none";
	public static PeptideModification phosphorylation=new PeptideModification("Phosphorylation (STY)", "Phosphorylation", 79.966331, new double[] {97.976896, 97.976896, 0.0}, new char[] {'S', 'T', 'Y'});
	public static PeptideModification acetylation=new PeptideModification("Acetylation (K)", "Acetylation", 42.010565, new double[1], new char[] {'K'});
	public static PeptideModification oxidation=new PeptideModification("Oxidation (MW)", "Oxidation", 15.994915, new double[2], new char[] {'M', 'W'});
	public static PeptideModification methylation=new PeptideModification("N-Methylation (KR)", "Methylation", 14.015650, new double[2], new char[] {'K', 'R'});
	public static PeptideModification ubiquitination=new PeptideModification("Ubiquitination (K)", "Ubiquitination", 114.042927, new double[1], new char[] {'K'});
	public static PeptideModification oglcnac=new PeptideModification("O-HexNAc (ST)", "OHexNAc", 203.079373, new double[] {203.079373, 203.079373}, new char[] {'S', 'T'});
	public static PeptideModification polymorphism=new PeptideModification("Nucleotide Polymorphism", "Polymorphism", 0.0, new double["ACDEFGHIKLMNPQRSTVWY".length()], "ACDEFGHIKLMNPQRSTVWY".toCharArray());

	public static final PeptideModification[] MODIFICATIONS=new PeptideModification[] {
			PeptideModification.phosphorylation,
			PeptideModification.acetylation,
			PeptideModification.oxidation,
			PeptideModification.methylation,
			PeptideModification.ubiquitination,
			PeptideModification.oglcnac
	};
	
	public static PeptideModification getModification(String name) {
		if (NO_MODIFICATION_NAME.equalsIgnoreCase(name)) return null;
		
		for (int i=0; i<MODIFICATIONS.length; i++) {
			if (MODIFICATIONS[i].getShortname().equalsIgnoreCase(name)) return MODIFICATIONS[i];
		}
		throw new EncyclopediaException("Sorry, only ["+getShortnameList()+"] are supported localization modifications.");
	}
	
	public static String getShortnameList() {
		StringBuilder sb=new StringBuilder();
		for (int i=0; i<MODIFICATIONS.length; i++) {
			if (i>0) sb.append(", ");
			sb.append(MODIFICATIONS[i].getShortname());
		}
		return sb.toString();
	}
	
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
	
	public static PeptideModification getModification(char aa, double modificationMass) {
		for (int i=0; i<MODIFICATIONS.length; i++) {
			if (MODIFICATIONS[i].isModificationMass(aa, modificationMass)) {
				return MODIFICATIONS[i];
			}
		}
		return null;
	}
	
	public static double getNeutralLoss(char aa, double modificationMass) {
		//if (true) return 0.0;
		PeptideModification mod=getModification(aa, modificationMass);
		if (mod==null) return 0.0;
		return mod.getNeutralLoss(aa);
	}
}
