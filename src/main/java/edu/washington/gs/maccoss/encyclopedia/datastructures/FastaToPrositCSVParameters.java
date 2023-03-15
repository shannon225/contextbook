package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;

import java.util.HashMap;

public class FastaToPrositCSVParameters {
	public static final int MAX_DEFAULT_NCE = 40;
	public static final int MIN_DEFAULT_NCE = 25;
	public static final int DEFAULT_DEFAULT_NCE = 33;

	public static final byte MAX_DEFAULT_CHARGE = 4;
	public static final byte MIN_DEFAULT_CHARGE = 1;
	public static final byte DEFAULT_DEFAULT_CHARGE = 3;

	public static final byte MIN_CHARGE = 1;
	public static final byte MAX_CHARGE = 6;
	public static final byte DEFAULT_MIN_CHARGE = 2;
	public static final byte DEFAULT_MAX_CHARGE = 3;

	public static final int MIN_MAX_MISSED_CLEAVAGE = 0;
	public static final int MAX_MAX_MISSED_CLEAVAGE = 3;
	public static final int DEFAULT_MAX_MISSED_CLEAVAGE = 1;

	public static final double MIN_MZ = 150.0;
	public static final double DEFAULT_MIN_MZ = 396.4;

	public static final double MAX_MZ = 1600.0;
	public static final double DEFAULT_MAX_MZ = 1002.7;
	
	public static final boolean DEFAULT_ADJUST_NCE_FOR_DIA = true;
	public static final boolean DEFAULT_ADD_DECOYS = false;

	public static final DigestionEnzyme DEFAULT_ENZYME = DigestionEnzyme.getEnzyme("Trypsin");

	protected final int defaultNCE;
	protected final byte defaultCharge;
	protected final byte minCharge;
	protected final byte maxCharge;
	protected final int maxMissedCleavage;
	protected final double minMz;
	protected final double maxMz;
	protected final DigestionEnzyme enzyme;
	protected final boolean adjustNCEForDIA;
	protected final boolean addDecoys;

	public FastaToPrositCSVParameters(int defaultNCE, byte defaultCharge, byte minCharge, byte maxCharge,
                                      int maxMissedCleavage, double minMz, double maxMz, DigestionEnzyme enzyme, boolean adjustNCEForDIA, boolean addDecoys) {
		this.defaultNCE = defaultNCE;
		this.defaultCharge = defaultCharge;
		this.minCharge = minCharge;
		this.maxCharge = maxCharge;
		this.maxMissedCleavage = maxMissedCleavage;
		this.minMz = minMz;
		this.maxMz = maxMz;
		this.enzyme = enzyme;
		this.adjustNCEForDIA=adjustNCEForDIA;
		this.addDecoys=addDecoys;
	}

	public String toString() {
		final StringBuilder sb=new StringBuilder();
		sb.append(" -defaultNCE ").append(defaultNCE).append("\n");
		sb.append(" -defaultCharge ").append(defaultCharge).append("\n");
		sb.append(" -minCharge ").append(minCharge).append("\n");
		sb.append(" -maxCharge ").append(maxCharge).append("\n");
		sb.append(" -maxMissedCleavage ").append(maxMissedCleavage).append("\n");
		sb.append(" -minMz ").append(minMz).append("\n");
		sb.append(" -maxMz ").append(maxMz).append("\n");
		sb.append(" -enzyme ").append(enzyme.getName()).append("\n");
		sb.append(" -adjustNCEForDIA ").append(adjustNCEForDIA).append("\n");
		sb.append(" -addDecoys ").append(addDecoys).append("\n");
		return sb.toString();
	}
	
	public HashMap<String, String> toParameterMap() {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("-defaultNCE", Integer.toString(defaultNCE));
		map.put("-defaultCharge", Byte.toString(defaultCharge));
		map.put("-minCharge", Integer.toString(minCharge));
		map.put("-maxCharge", Integer.toString(maxCharge));
		map.put("-maxMissedCleavage", Integer.toString(maxMissedCleavage));
		map.put("-minMz", Double.toString(minMz));
		map.put("-maxMz", Double.toString(maxMz));
		map.put("-enzyme", enzyme.getName());
		map.put("-adjustNCEForDIA", Boolean.toString(adjustNCEForDIA));
		map.put("-addDecoys", Boolean.toString(addDecoys));
		return map;
	}

	public int getDefaultNCE() {
		return defaultNCE;
	}

	public byte getDefaultCharge() {
		return defaultCharge;
	}

	public byte getMinCharge() {
		return minCharge;
	}

	public byte getMaxCharge() {
		return maxCharge;
	}

	public int getMaxMissedCleavage() {
		return maxMissedCleavage;
	}

	public double getMinMz() {
		return minMz;
	}

	public double getMaxMz() {
		return maxMz;
	}

	public DigestionEnzyme getEnzyme() {
		return enzyme;
	}
	
	public boolean isAddDecoys() {
		return addDecoys;
	}
	
	public boolean isAdjustNCEForDIA() {
		return adjustNCEForDIA;
	}
}