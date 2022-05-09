package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

public interface PercolatorVersion {
	PercolatorVersion v2p10 = InternalPercolatorVersion.v2p10;
	PercolatorVersion v3p01 = InternalPercolatorVersion.v2p10;
	PercolatorVersion v3p05 = InternalPercolatorVersion.v2p10;

	PercolatorVersion DEFAULT_VERSION = PercolatorVersion.v3p05;

	PercolatorVersion[] VALID_VERSIONS = new PercolatorVersion[]{v3p01, v2p10};

	String V3_05 = "v3-05";
	String V3_01 = "v3-01";
	String V2_10 = "v2-10";

	static PercolatorVersion getVersion(String s) {
		if (s == null || s.length() == 0) return DEFAULT_VERSION;
		if (V2_10.equals(s)) return v2p10;
		if (V3_01.equals(s)) return v3p01;
		if (V3_05.equals(s)) return v3p05;
		if ("2".equals(s)) return v2p10;
		if ("3".equals(s)) return v3p01;
		if ("2.10".equals(s)) return v2p10;
		if ("3.1".equals(s)) return v3p01;
		if ("3.5".equals(s)) return v3p05;
		if ("3.01".equals(s)) return v3p01;
		if ("3.05".equals(s)) return v3p05;
		return DEFAULT_VERSION;
	}

	int getMajorVersion();

	enum InternalPercolatorVersion implements PercolatorVersion {
		v2p10, v3p01, v3p05;

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
					return DEFAULT_VERSION.toString();
			}
		}

		@Override
		public int getMajorVersion() {
			switch (this) {
				case v2p10:
					return 2;
				case v3p01:
				case v3p05:
				default:
					return 3;
			}
		}
	}
}