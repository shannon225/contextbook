package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;

public class RetentionTimeBoundary extends Range implements HasRetentionTime {
	private final float apexRT;

	public RetentionTimeBoundary(float minRT, float apexRT, float maxRT) {
		super(minRT, maxRT);
		this.apexRT = apexRT;
	}

	public static RetentionTimeBoundary getMedianBoundaries(float[] rts) {
		float[] clone=rts.clone();
		Arrays.sort(clone);
		return new RetentionTimeBoundary(clone[0], QuickMedian.median(clone), clone[clone.length-1]);
	}
	
	public float getMinRT() {
		return getStart();
	}

	public float getMaxRT() {
		return getStop();
	}
	
	public float getApexRT() {
		return apexRT;
	}

	@Override
	public float getRetentionTimeInSec() {
		return apexRT;
	}
	
	public boolean contains(float rt) {
		return (rt>=getMinRT()&&rt<=getMaxRT());
	}
	
	@Override
	public String toString() {
		return "["+getMinRT()+" - "+apexRT+" - "+getMaxRT()+"]";
	}
}
