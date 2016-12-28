package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.Comparator;

public class PSMDataMassComparator implements Comparator<PSMData> {
	public PSMDataMassComparator() {
	}

	@Override
	public int compare(PSMData o1, PSMData o2) {
		if (o1==null&&o2==null) {
			return 0;
		} else if (o1==null) {
			return -1;
		} else if (o2==null) {
			return 1;
		}
		
		return Double.compare(o1.getUnchargedMass(), o2.getUnchargedMass());
	}
}
