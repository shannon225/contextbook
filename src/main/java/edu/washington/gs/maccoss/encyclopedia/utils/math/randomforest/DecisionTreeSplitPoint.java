package edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest;

import java.util.ArrayList;

public class DecisionTreeSplitPoint implements Comparable<DecisionTreeSplitPoint> {
	public final int index;
	public final float cutoff;
	public final float gain;

	public DecisionTreeSplitPoint(int index, float cutoff, float gain) {
		this.index=index;
		this.cutoff=cutoff;
		this.gain=gain;
	}
	
	@Override
	public int compareTo(DecisionTreeSplitPoint o) {
		if (o==null) return 1;

		int c=Integer.compare(index, o.index);
		if (c!=0) return 0;
		c=Float.compare(cutoff, o.cutoff);
		if (c!=0) return 0;
		c=Float.compare(gain, o.gain);
		if (c!=0) return 0;
		return 0;
	}

	public boolean isRight(Record r) {
		if (r.getData()[index]>cutoff) {
			return true;
		} else {
			return false;
		}
	}

	@Override
	public String toString() {
		return index+","+cutoff;
	}

	@SuppressWarnings("unchecked")
	public ArrayList<Record>[] split(ArrayList<Record> records) {
		ArrayList<Record> left=new ArrayList<Record>();
		ArrayList<Record> right=new ArrayList<Record>();
		for (Record record : records) {
			if (isRight(record)) {
				right.add(record);
			} else {
				left.add(record);
			}
		}
		return new ArrayList[] { left, right };
	}
}
