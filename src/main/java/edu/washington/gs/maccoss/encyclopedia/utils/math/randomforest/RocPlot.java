package edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest;

import gnu.trove.list.array.TFloatArrayList;

public class RocPlot {
	TFloatArrayList fprs=new TFloatArrayList();
	TFloatArrayList tprs=new TFloatArrayList();

	public RocPlot() {
		fprs.add(0.0f);
		tprs.add(0.0f);
	}

	/**
	 * assumes adding in order of fpr
	 * 
	 * @param fpr
	 * @param tpr
	 */
	public void addData(float fpr, float tpr) {
		if (fprs.size()>0&&fpr==fprs.get(fprs.size()-1)) {
			// overwrite lower tpr
			tprs.set(tprs.size()-1, Math.max(tprs.get(tprs.size()-1), tpr));
		} else {
			// increment
			fprs.add(fpr);
			tprs.add(tpr);
		}
	}

	public String toString() {
		StringBuilder sb=new StringBuilder();
		for (int i=0; i<fprs.size(); i++) {
			sb.append(fprs.get(i)+"\t"+tprs.get(i)+"\n");
		}
		return sb.toString();
	}

	public float getAUC() {
		float sum=0.0f;
		for (int i=0; i<fprs.size(); i++) {
			if (i==fprs.size()-1) {
				sum+=getTrapezoidArea(fprs.get(i), tprs.get(i), 1.0f, 1.0f);
			} else {
				sum+=getTrapezoidArea(fprs.get(i), tprs.get(i), fprs.get(i+1), tprs.get(i+1));
			}
		}
		return sum;
	}

	public float getTPR(float fpr) {
		if (fpr>=1.0f)
			return 1.0f;
		int index=fprs.binarySearch(fpr);
		index=index>=0 ? index : -(index+1);
		if (index<1)
			return tprs.get(0);
		if (index>=tprs.size())
			return tprs.get(tprs.size()-1);

		return interpolate(fprs.get(index-1), fpr, fprs.get(index), tprs.get(index-1), tprs.get(index));
	}

	private float interpolate(float minX, float x, float maxX, float minY, float maxY) {
		return ((maxY-minY)/((maxX-minX)))*((x-minX))+minY;
	}

	private float getTrapezoidArea(float x1, float y1, float x2, float y2) {
		float deltaX=x2-x1;
		float square=deltaX*y1;
		float triangle=deltaX*(y2-y1)/2.0f;
		return square+triangle;
	}
}
