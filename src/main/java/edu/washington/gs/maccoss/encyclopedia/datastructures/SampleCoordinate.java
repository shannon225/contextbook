package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.algorithms.EncyclopediaElibParser.Coordinate;

public class SampleCoordinate {
	private final int replicate;
	private final int sample;

	public SampleCoordinate(int replicate, int sample) {
		this.replicate=replicate;
		this.sample=sample;
	}

	public int getReplicate() {
		return replicate;
	}

	public int getSample() {
		return sample;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj==null) return false;
		if (!(obj instanceof Coordinate)) return false;
		return hashCode()==obj.hashCode();
	}

	@Override
	public int hashCode() {
		return sample+replicate*10000;
	}
}
