package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

public class SampleCoordinate {
	private final int replicateIndex;
	private final int sampleIndex;

	public SampleCoordinate(int replicateIndex, int sampleIndex) {
		this.replicateIndex=replicateIndex;
		this.sampleIndex=sampleIndex;
	}

	public int getReplicateIndex() {
		return replicateIndex;
	}

	public int getSampleIndex() {
		return sampleIndex;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj==null) return false;
		if (!(obj instanceof SampleCoordinate)) return false;
		if (hashCode()==obj.hashCode()) {
			SampleCoordinate c=(SampleCoordinate)obj;
			return (replicateIndex==c.replicateIndex&&sampleIndex==c.sampleIndex);
		} else {
			return false;
		}
	}

	@Override
	public int hashCode() {
		return sampleIndex+replicateIndex*10000;
	}
}
