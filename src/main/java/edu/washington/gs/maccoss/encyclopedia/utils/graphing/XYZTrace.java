package edu.washington.gs.maccoss.encyclopedia.utils.graphing;

import java.util.ArrayList;

import org.jfree.data.DomainOrder;
import org.jfree.data.general.DatasetChangeListener;
import org.jfree.data.general.DatasetGroup;
import org.jfree.data.xy.XYZDataset;

public class XYZTrace implements XYZDataset {
	private final String name;
	private final ArrayList<XYZPoint> points;
	private final float maxZ;
	private final float minZ;

	public XYZTrace(String name, ArrayList<XYZPoint> points) {
		this.name=name;
		this.points=points;
		
		double localMaxZ=-Float.MAX_VALUE;
		double localMinZ=Float.MAX_VALUE;
		for (XYZPoint point : points) {
			if (point.z>localMaxZ) localMaxZ=point.z;
			if (point.z<localMinZ) localMinZ=point.z;
		}
		maxZ=(float)localMaxZ;
		minZ=(float)localMinZ;
	}
	
	public float getMaxZ() {
		return maxZ;
	}
	public float getMinZ() {
		return minZ;
	}

	@Override
	public Number getX(int series, int item) {
		return points.get(item).getX();
	}

	@Override
	public double getXValue(int series, int item) {
		return points.get(item).getX();
	}

	@Override
	public Number getY(int series, int item) {
		return points.get(item).getY();
	}

	@Override
	public double getYValue(int series, int item) {
		return points.get(item).getY();
	}

	@Override
	public Number getZ(int series, int item) {
		return points.get(item).getZ();
	}

	@Override
	public double getZValue(int series, int item) {
		return points.get(item).getZ();
	}

	@Override
	public int getSeriesCount() {
		return 1;
	}

	@Override
	public int getItemCount(int series) {
		return points.size();
	}

	@Override
	public void addChangeListener(DatasetChangeListener listener) {
		// ignore - this dataset never changes
	}

	@Override
	public void removeChangeListener(DatasetChangeListener listener) {
		// ignore
	}

	@Override
	public DatasetGroup getGroup() {
		return null;
	}

	@Override
	public void setGroup(DatasetGroup group) {
		// ignore
	}

	@SuppressWarnings("rawtypes")
	@Override
	public Comparable getSeriesKey(int series) {
		return name;
	}

	@SuppressWarnings("rawtypes")
	@Override
	public int indexOf(Comparable seriesKey) {
		return 0;
	}

	@Override
	public DomainOrder getDomainOrder() {
		return DomainOrder.ASCENDING;
	}

}
