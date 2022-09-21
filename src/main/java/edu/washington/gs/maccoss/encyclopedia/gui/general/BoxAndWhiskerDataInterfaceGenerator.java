package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.util.List;

import org.jfree.data.statistics.BoxAndWhiskerCategoryDataset;
import org.jfree.data.statistics.BoxAndWhiskerXYDataset;

public class BoxAndWhiskerDataInterfaceGenerator {
	public static BoxAndWhiskerDataInterface getInterface(final BoxAndWhiskerXYDataset dataset) {
		return new BoxAndWhiskerDataInterface() {
			@Override
			public Number getQ3Value(int row, int column) {
				return dataset.getQ3Value(row, column);
			}
			
			@Override
			public Number getQ1Value(int row, int column) {
				return dataset.getQ1Value(row, column);
			}
			
			@Override
			public List getOutliers(int row, int column) {
				return dataset.getOutliers(row, column);
			}
			
			@Override
			public Number getMinRegularValue(int row, int column) {
				return dataset.getMinRegularValue(row, column);
			}
			
			@Override
			public Number getMedianValue(int row, int column) {
				return dataset.getMedianValue(row, column);
			}
			
			@Override
			public Number getMaxRegularValue(int row, int column) {
				return dataset.getMaxRegularValue(row, column);
			}
		};
	}
	
	public static BoxAndWhiskerDataInterface getInterface(final BoxAndWhiskerCategoryDataset dataset) {
		return new BoxAndWhiskerDataInterface() {
			@Override
			public Number getQ3Value(int row, int column) {
				return dataset.getQ3Value(row, column);
			}
			
			@Override
			public Number getQ1Value(int row, int column) {
				return dataset.getQ1Value(row, column);
			}
			
			@Override
			public List getOutliers(int row, int column) {
				return dataset.getOutliers(row, column);
			}
			
			@Override
			public Number getMinRegularValue(int row, int column) {
				return dataset.getMinRegularValue(row, column);
			}
			
			@Override
			public Number getMedianValue(int row, int column) {
				return dataset.getMedianValue(row, column);
			}
			
			@Override
			public Number getMaxRegularValue(int row, int column) {
				return dataset.getMaxRegularValue(row, column);
			}
		};
	}

}
