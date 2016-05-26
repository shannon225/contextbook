package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.util.StringTokenizer;

public class Version implements Comparable<Version> {
	private final int major;
	private final int minor;
	private final int revision;

	public Version(int major, int minor, int revision) {
		this.major=major;
		this.minor=minor;
		this.revision=revision;
	}

	public Version(String versionString) {
		if (versionString==null) {
			major=0;
			minor=0;
			revision=0;
		} else {
			StringTokenizer st=new StringTokenizer(versionString, ".");
			major=Integer.parseInt(st.nextToken());
			minor=Integer.parseInt(st.nextToken());
			if (st.hasMoreTokens()) {
				revision=Integer.parseInt(st.nextToken());
			} else {
				revision=0;
			}
		}
	}

	public String toString() {
		if (revision==0) {
			return major+"."+minor;
		} else {
			return major+"."+minor+"."+revision;
		}
	}
	
	public boolean amIAbove(Version v) {
		return compareTo(v)>0;
	}

	@Override
	public int compareTo(Version o) {
		if (o==null) return 1;
		int c=Integer.compare(major, o.major);
		if (c!=0) return c;
		c=Integer.compare(minor, o.minor);
		if (c!=0) return c;
		c=Integer.compare(revision, o.revision);
		return c;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof Version) {
			return compareTo((Version)obj)==0;
		} else {
			return false;
		}
	}

	@Override
	public int hashCode() {
		return major+minor*1000+revision*1000000;
	}
}
