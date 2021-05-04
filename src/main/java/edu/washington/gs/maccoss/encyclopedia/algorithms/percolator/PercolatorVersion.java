package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

public enum PercolatorVersion {
	v2p10,v3p01, v3p05;
	
	public static final String V3_05="v3-05";
	public static final String V3_01="v3-01";
	public static final String V2_10="v2-10";
	public static final PercolatorVersion DEFAULT_VERSION=PercolatorVersion.v3p01;
	public static final PercolatorVersion[] VALID_VERSIONS=new PercolatorVersion[] {v3p05, v3p01, v2p10};

	public static PercolatorVersion getVersion(String s) {
		if (s==null||s.length()==0) return v3p05;
		if (V2_10.equals(s)) return v2p10;
		if (V3_01.equals(s)) return v3p01;
		if (V3_05.equals(s)) return v3p05;
		if ("2".equals(s)) return v2p10;
		if ("3".equals(s)) return v3p05;
		if ("2.10".equals(s)) return v2p10;
		if ("3.1".equals(s)) return v3p01;
		if ("3.5".equals(s)) return v3p05;
		if ("3.01".equals(s)) return v3p01;
		if ("3.05".equals(s)) return v3p05;
		return v3p05;
	}
	
	@Override
	public String toString() {
		switch (this) {
		case v2p10:
			return V2_10;
		case v3p01:
			return V3_01;
		case v3p05:
			return V3_05;
		default:
			return V3_05;
		}
	}
	
	public int getMajorVersion() {
		switch (this) {
		case v2p10:
			return 2;
		case v3p01:
			return 3;
		case v3p05:
			return 3;
		default:
			return 3;
		}
	}
}
