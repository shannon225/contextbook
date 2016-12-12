package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.StringTokenizer;

public class FragmentIon implements Comparable<FragmentIon> {
	private static final String ARCHIVE_DELIMINATOR=";";
	private static final MassTolerance tolerance=new MassTolerance(0.1f); 
	public final double mass;
	public final byte index;
	public final IonType type;

	public FragmentIon(double mass, byte index, IonType type) {
		this.mass=mass;
		this.index=index;
		this.type=type;
	}
	
	public FragmentIon increment(double deltaMass) {
		return new FragmentIon(mass+deltaMass, (byte)(index+1), type);
	}
	
	public FragmentIon neutralLoss(double deltaMass) {
		return new FragmentIon(mass-deltaMass, index, IonType.getNL(type));
	}
	
	public IonType getType() {
		return type;
	}
	
	public String toCanonicalIonTypeString() {
		return IonType.toString(IonType.getCanonicalIonType(type), index);
	}
	
	@Override
	public String toString() {
		return IonType.toString(type, index);
	}
	
	public String toArchiveString() {
		return IonType.toString(type)+ARCHIVE_DELIMINATOR+index;
	}
	
	public FragmentIon fromArchiveString(String s, double mass) {
		StringTokenizer st=new StringTokenizer(s, ARCHIVE_DELIMINATOR);
		IonType type=IonType.fromString(st.nextToken());
		byte index=Byte.parseByte(st.nextToken());
		return new FragmentIon(mass, index, type);
	}
	
	@Override
	public int compareTo(FragmentIon o) {
		if (o==null) return 1;
		return tolerance.compareTo(mass, o.mass);
	}
	
	@Override
	public int hashCode() {
		return Double.hashCode(mass);
	}
	
	@Override
	public boolean equals(Object o) {
		if (o==null) return false;
		return tolerance.compareTo(mass, ((FragmentIon)o).mass)==0;
	}
	
	public static double[] getMasses(FragmentIon[] ions) {
		double[] masses=new double[ions.length];
		for (int i=0; i<masses.length; i++) {
			masses[i]=ions[i].mass;
		}
		return masses;
	}
}
