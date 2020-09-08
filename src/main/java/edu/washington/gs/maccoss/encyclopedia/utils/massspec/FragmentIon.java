package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;

import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;

public final class FragmentIon implements Comparable<FragmentIon> {
	private static final String INDEX_DELIMITER = ";";
	private static final String ARCHIVE_DELIMITER = "|";

	private final double mass;
	private final byte index;
	private final IonType type;

	public FragmentIon(double mass, byte index, IonType type) {
		this.mass = mass;
		this.index = index;
		this.type = type;
	}
	
	public double getMass() {
		return mass;
	}
	
	public byte getIndex() {
		return index;
	}

	public IonType getType() {
		return type;
	}

	public String getName() {
		return IonType.toString(type, index);
	}

	public static String toArchiveString(FragmentIon[] ions) {
		StringBuilder sb = new StringBuilder();
		for (FragmentIon ion : ions) {
			if (sb.length() > 0) {
				sb.append(ARCHIVE_DELIMITER);
			}
			sb.append(IonType.toString(ion.type));
			sb.append(INDEX_DELIMITER);
			sb.append(ion.index);
			sb.append(INDEX_DELIMITER);
			sb.append(ion.mass);
		}
		return sb.toString();
	}

	public static FragmentIon[] fromArchiveString(String s) {
		if (s == null || s.trim().length() == 0 || s.equalsIgnoreCase("null")) {
			return new FragmentIon[0];
		}
		StringTokenizer st = new StringTokenizer(s, ARCHIVE_DELIMITER);
		ArrayList<FragmentIon> ions = new ArrayList<>();
		while (st.hasMoreTokens()) {
			StringTokenizer st2 = new StringTokenizer(st.nextToken(), INDEX_DELIMITER);
			String tag=st2.nextToken();
			IonType type = IonType.fromString(tag);
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

	//TODO: write test cases for this functionality
	public static FragmentIon[] getUniqueFragments(FragmentIon[] fragments, MassTolerance tolerance) {
		final List<FragmentIon> work = Lists.newArrayList();

		Arrays.stream(fragments)
				.sorted(Ordering.natural().onResultOf(FragmentIon::getType)) // sort by declaration order of fragment types -- this is roughly by priority
				.forEach(ion -> {
					boolean exists = false;
					for (FragmentIon existing : work) {
						if (tolerance.equals(ion.mass, existing.mass)) {
							exists = true;
							break;
						}
					}

					if (!exists) {
						work.add(ion);
					}
				});

		return work.stream()
				.sorted() // sort by natural ordering of ions
				.toArray(FragmentIon[]::new);
	}

	public FragmentIon increment(double deltaMass) {
		return new FragmentIon(mass + deltaMass, (byte) (index + 1), type);
	}

	public FragmentIon neutralLoss(double deltaMass) {
		return new FragmentIon(mass - deltaMass, index, IonType.getNL(type));
	}

	public String toCanonicalIonTypeString() {
		return IonType.toString(IonType.getCanonicalIonType(type), index);
	}

	@Override
	public int hashCode() {
		// Note that equal objects will always have identical masses (see below)
		return (int)(mass*100.0);
	}

	@Override
	public boolean equals(Object o) {
		return o instanceof FragmentIon && compareTo((FragmentIon) o) == 0;
	}

	@Override
	public String toString() {
		return getName();
	}
	
	public Color getColor() {
		return RandomGenerator.randomColor(toString().hashCode());
	}

	@Override
	public int compareTo(FragmentIon o) {
		if (o == null) {
			return 1;
		}
		int c=tolerance.compareTo(mass, o.mass);
		if (c!=0) return c;

		// Comparison uses exact mass as well as type and index. Natural ordering will be by mass, with ties settled
		// by type (ordered by declaration order), then index. This will be transitive and consistent with equals, as
		// all comparisons are exact.
		return ComparisonChain.start()
				.compare(type, o.type)
				.compare(index, o.index)
				.result();
	}
	
	private static final MassTolerance tolerance=new MassTolerance(0.1); // 1 ppm is about the accuracy of floats 
}
