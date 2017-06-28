package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.StringTokenizer;

public class FragmentIon implements Comparable<FragmentIon> {
	private static final String INDEX_DELIMINATOR=";";
	private static final String ARCHIVE_DELIMINATOR="|";
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
	
	public static String toArchiveString(FragmentIon[] ions) {
		StringBuilder sb=new StringBuilder();
		for (FragmentIon ion : ions) {
			if (sb.length()>0) {
				sb.append(ARCHIVE_DELIMINATOR);
			}
			sb.append(IonType.toString(ion.type));
			sb.append(INDEX_DELIMINATOR);
			sb.append(ion.index);
			sb.append(INDEX_DELIMINATOR);
			sb.append(ion.mass);
		}
		return sb.toString();
	}
	
	public static FragmentIon[] fromArchiveString(String s) {
		if (s==null||s.trim().length()==0) {
			return new FragmentIon[0];
		}
		StringTokenizer st=new StringTokenizer(s, ARCHIVE_DELIMINATOR);
		ArrayList<FragmentIon> ions=new ArrayList<FragmentIon>();
		while (st.hasMoreTokens()) {
			StringTokenizer st2=new StringTokenizer(st.nextToken(), INDEX_DELIMINATOR);
			IonType type=IonType.fromString(st2.nextToken());
			byte index=Byte.parseByte(st2.nextToken());
			double mass=Double.parseDouble(st2.nextToken());
			ions.add(new FragmentIon(mass, index, type));
		}
		return ions.toArray(new FragmentIon[ions.size()]);
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
