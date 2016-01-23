package edu.washington.gs.maccoss.encyclopedia.utils.math;

public interface Function {
	public float getYValue(float xi);
	public boolean isXInsideBoundaries(float xi);
	public boolean isYInsideBoundaries(float yi);
	public float getXValue(float yi);
}
