package edu.washington.gs.maccoss.encyclopedia.utils.math;

public class QuickMedianInteger {
	public static void main(String[] args) {
		int[] data=new int[] {368149,219199,68060,105270};

		int median=QuickMedianInteger.select(data, 0.5f);
		int q1=QuickMedianInteger.select(data, 0.25f);
		int q3=QuickMedianInteger.select(data, 0.75f);
		int lowerOutlierThreshold=QuickMedianInteger.select(data, 0.05f);
		int upperOutlierThreshold=QuickMedianInteger.select(data, 0.95f);

		System.out.println(General.toString(data));
		System.out.println(General.toString(new int[] {lowerOutlierThreshold, q1, median, q3, upperOutlierThreshold}));
	}
	
	public static int median(int[] data) {
		return select(data, 0.5f);
	}
	
	public static int iqr(int[] data) {
		return select(data, 0.75f)-select(data, 0.25f);
	}
	
	public static int range90(int[] data) {
		return select(data, 0.95f)-select(data, 0.05f);
	}

	/**
	 * this is destructive to the array order! (it sorts it)
	 * @param data
	 * @param desiredPercentage
	 * @return
	 */
	public static int select(int[] data, float desiredPercentage) {
		if (data.length==0) return 0;
		if (data.length==1) return data[0];

		// The exact position of the value we want to select. Note that this might fall between two values.
		final float pos = (data.length - 1) * desiredPercentage;

		// The fractional part of the position gives the proportion we should use to combine values if we fall between two.
		final float prop = pos % 1;

		// The desired target -- one more than the integral part of the position.
		int targetIndex=Math.round(pos - prop) + 1;

		int left=0;
		int right=data.length-1;
		
		int seed=RandomGenerator.randomInt(data.length);

		while (true) {
			seed=RandomGenerator.randomInt(seed);
			if (left==right) {
				if (left < data.length - 1) {
					// guaranteed to be locally sorted because this point will always be sandwiched between two pivot points
					return Math.round(((1 - prop) * data[left]) + (prop * data[left + 1]));
				} else {
					return data[left];
				}
			}
			int pivotIndex=left+(int)Math.round((right-left)*RandomGenerator.floatFromRandomInt(seed));
			if (pivotIndex<0) {
				pivotIndex=0;
			}
			int pivotNewIndex=partition(data, left, right, pivotIndex);

			int pivotDist=pivotNewIndex-left+1;
			if (pivotDist==targetIndex) {
				if (pivotNewIndex<(data.length-1)) {
					float lowest=data[pivotNewIndex+1];
					for (int i=pivotNewIndex+2; i<=right; i++) {
						if (data[i]<lowest) {
							lowest=data[i];
						}
					}
					return Math.round(((1 - prop) * data[pivotNewIndex]) + (prop * lowest));
				} else {
					return data[pivotNewIndex];
				}
			} else if (targetIndex<pivotDist) {
				// the target is left of the pivot
				if (pivotNewIndex < 1) {
					// There's no data left of the pivot -- this means left is zero (because pivotNewIndex >= left always).
					// We just set right to zero and let the loop terminate next time around (because left==right).
					right = 0;
				} else {
					right = pivotNewIndex - 1;
				}
			} else {
				// the target is right of the pivot
				targetIndex=targetIndex-pivotDist;
				left=pivotNewIndex+1;
			}
		}
	}

	private static int partition(int[] data, int left, int right, int pivotIndex) {
		int pivotValue=data[pivotIndex];
		swap(data, pivotIndex, right);
		int storeIndex=left;
		for (int i=left; i<=right; i++) {
			if (data[i]<pivotValue) {
				swap(data, storeIndex, i);
				storeIndex++;
			}
		}
		swap(data, right, storeIndex);
		return storeIndex;
	}

	/**
	 * Method to swap to elements in an array.
	 * 
	 * @param a
	 *            an array of objects.
	 * @param index1
	 *            the index of the first object.
	 * @param index2
	 *            the index of the second object.
	 */
	public static final void swap(int[] a, int index1, int index2) {
		int tmp=a[index1];
		a[index1]=a[index2];
		a[index2]=tmp;
	}
}
