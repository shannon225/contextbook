package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.util.List;

public interface BoxAndWhiskerDataInterface {
	public Number getQ1Value(int row, int column);
	public Number getQ3Value(int row, int column);
	public Number getMaxRegularValue(int row, int column);
	public Number getMinRegularValue(int row, int column);
	public Number getMedianValue(int row, int column);
	public List getOutliers(int row, int column);
}
