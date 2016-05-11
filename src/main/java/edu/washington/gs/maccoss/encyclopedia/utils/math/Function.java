package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;

public interface Function {
	public float getYValue(float xi);
	public boolean isXInsideBoundaries(float xi);
	public boolean isYInsideBoundaries(float yi);
	public float getXValue(float yi);
	public ArrayList<XYPoint> getKnots();
}
