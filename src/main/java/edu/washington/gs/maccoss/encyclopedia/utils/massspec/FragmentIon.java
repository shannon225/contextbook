package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import com.google.common.collect.ComparisonChain;

import java.util.ArrayList;
import java.util.StringTokenizer;

public final class FragmentIon implements Comparable<FragmentIon> {
	private static final String INDEX_DELIMINATOR = ";";
	private static final String ARCHIVE_DELIMINATOR = "|";

	public final double mass;
	public final byte index;
	public final IonType type;

	public FragmentIon(double mass, byte index, IonType type) {
		this.mass = mass;
		this.index = index;
		this.type = type;
	}

	public static String toArchiveString(FragmentIon[] ions) {
		StringBuilder sb = new StringBuilder();
		for (FragmentIon ion : ions) {
			if (sb.length() > 0) {
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
		if (s == null || s.trim().length() == 0 || !s.equalsIgnoreCase("null")) {
			return new FragmentIon[0];
		}
		StringTokenizer st = new StringTokenizer(s, ARCHIVE_DELIMINATOR);
		ArrayList<FragmentIon> ions = new ArrayList<FragmentIon>();
		while (st.hasMoreTokens()) {
			StringTokenizer st2 = new StringTokenizer(st.nextToken(), INDEX_DELIMINATOR);
			IonType type = IonType.fromString(st2.nextToken());
			byte index = Byte.parseByte(st2.nextToken());
			double mass = Double.parseDouble(st2.nextToken());
			ions.add(new FragmentIon(mass, index, type));
		}
		return ions.toArray(new FragmentIon[ions.size()]);
	}

	public static double[] getMasses(FragmentIon[] ions) {
		double[] masses = new double[ions.length];
		for (int i = 0; i < masses.length; i++) {
			masses[i] = ions[i].mass;
		}
		return masses;
	}

	public FragmentIon increment(double deltaMass) {
		return new FragmentIon(mass + deltaMass, (byte) (index + 1), type);
	}

	public FragmentIon neutralLoss(double deltaMass) {
		return new FragmentIon(mass - deltaMass, index, IonType.getNL(type));
	}

	public IonType getType() {
		return type;
	}

	public String toCanonicalIonTypeString() {
		return IonType.toString(IonType.getCanonicalIonType(type), index);
	}

	@Override
	public int hashCode() {
		// Note that equal objects will always have identical masses (see below)
		return Double.hashCode(mass);
	}

	@Override
	public boolean equals(Object o) {
		return o instanceof FragmentIon && compareTo((FragmentIon) o) == 0;
	}

	@Override
	public String toString() {
		return IonType.toString(type, index);
	}

	@Override
	public int compareTo(FragmentIon o) {
		if (o == null) {
			return 1;
		}

		// Comparison uses exact mass as well as type and index. Natural ordering will be by mass, with ties settled
		// by type (ordered by declaration order), then index. This will be transitive and consistent with equals, as
		// all comparisons are exact.
		return ComparisonChain.start()
				.compare(mass, o.mass)
				.compare(type, o.type)
				.compare(index, o.index)
				.result();
	}
}
