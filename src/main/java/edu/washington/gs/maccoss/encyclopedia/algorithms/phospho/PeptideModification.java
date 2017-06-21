package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class PeptideModification {
	public static final String NO_MODIFICATION_NAME="none";
	public static PeptideModification phosphorylation=new PeptideModification("Phosphorylation (STY)", "Phosphorylation", 79.966331, new char[] {'S', 'T', 'Y'});
	public static PeptideModification acetylation=new PeptideModification("Acetylation (K)", "Acetylation", 42.010565, new char[] {'K'});
	public static PeptideModification oxidation=new PeptideModification("Oxidation (MW)", "Oxidation", 15.994915, new char[] {'M', 'W'});
	public static PeptideModification methylation=new PeptideModification("N-Methylation (KR)", "Methylation", 14.015650, new char[] {'K', 'R'});
	public static PeptideModification ubiquitination=new PeptideModification("Ubiquitination (K)", "Ubiquitination", 114.042927, new char[] {'K'});
	public static PeptideModification oglcnac=new PeptideModification("O-GlcNAc (ST)", "OGlcNAc", 203.079373, new char[] {'S', 'T'});

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
	private final int nominalMass;
	private final char[] modifiableAAs;

	public PeptideModification(String name, String shortname, double mass, char[] modifiableAAs) {
		this.name=name;
		this.shortname=shortname;
		this.mass=mass;
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
}
