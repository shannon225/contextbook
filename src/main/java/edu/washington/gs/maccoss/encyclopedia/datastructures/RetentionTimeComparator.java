package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.Comparator;

public class RetentionTimeComparator implements Comparator<HasRetentionTime> {
	public static final RetentionTimeComparator comparator=new RetentionTimeComparator();
	
	@Override
	public int compare(HasRetentionTime o1, HasRetentionTime o2) {
		if (o1==null&&o2==null) return 0;
		else if (o1==null) return 1;
		else if (o2==null) return -1;
		
		return Float.compare(o1.getRetentionTimeInSec(), o2.getRetentionTimeInSec());
	}

	public static HasRetentionTime getComparable(final float retentionTimeInSec) {
		return new HasRetentionTime() {
			@Override
			public float getRetentionTimeInSec() {
				return retentionTimeInSec;
			}
		};
	}
}
