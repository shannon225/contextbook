package edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import gnu.trove.list.array.TFloatArrayList;

public class RocPlot {
	private static final int DEFAULT_ORDER = 4;
	private final float rounder;
	private final String name;
	private final TFloatArrayList fprs=new TFloatArrayList();
	private final TFloatArrayList tprs=new TFloatArrayList();


	public RocPlot() {
		this(DEFAULT_ORDER);
	}
	public RocPlot(String name) {
		this(DEFAULT_ORDER, name);
	}
	public RocPlot(int order) {
		this(order, "ROC Plot");
	}
	public RocPlot(int order, String name) {
		this.rounder=(float)Math.pow(10, order);
		this.name=name;
		fprs.add(0.0f);
		tprs.add(0.0f);
	}
	
	public int size() {
		return fprs.size();
	}

	/**
	 * assumes adding in order of fpr
	 * 
	 * @param fpr
	 * @param tpr
	 */
	public void addData(float fpr, float tpr) {
		fpr=Math.round(fpr*rounder)/rounder;
		
		if (fprs.size()>0&&fpr<=fprs.get(fprs.size()-1)) {
			// overwrite lower tpr
			fprs.set(fprs.size()-1, fpr);
			tprs.set(tprs.size()-1, Math.max(tprs.get(tprs.size()-1), tpr));

			while (true) {
				if (fprs.size()>1&&fpr<=fprs.get(fprs.size()-2)) {
					fprs.removeAt(fprs.size()-2);
					tprs.removeAt(tprs.size()-2);
				} else {
					break;
				}
			}
			
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
	
	public XYTrace getTrace() {
		return new XYTrace(fprs.toArray(), tprs.toArray(), GraphType.line, name);
	}

	public float getAUC() {
		return getAUC(1.0f);
	}
	public float getAUC(float maxFPR) {
		float sum=0.0f;
		for (int i=0; i<fprs.size(); i++) {
			if (fprs.get(i)>maxFPR) break;
			
			if (i==fprs.size()-1) {
				sum+=getTrapezoidArea(fprs.get(i), tprs.get(i), 1.0f, 1.0f);
			} else {
				sum+=getTrapezoidArea(fprs.get(i), tprs.get(i), fprs.get(i+1), tprs.get(i+1));
			}
		}
		return sum/maxFPR;
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
