package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.WeightedValue;
import junit.framework.TestCase;

public class KDETest extends TestCase {
	public static void main(String[] args) {
		KDE dist=new KDE(getData(), 1);
		Range range=new Range(-2, 8);
		Charter.launchChart(Charter.getChart(dist, range), "kde");
	}
	
	public void testKDE() {
		KDE dist=new KDE(getData(), 1);
		assertEquals(4, dist.getMode(), 0.1);
		assertEquals(3.7, dist.getMean(), 0.1);
		assertEquals(1.3, dist.getStdev(), 0.1);
	}

	public static ArrayList<WeightedValue> getData() {
		ArrayList<WeightedValue> data=new ArrayList<WeightedValue>();
		data.add(new WeightedValue(1, 0.5));
		data.add(new WeightedValue(2, 1));
		data.add(new WeightedValue(3, 0.5));
		data.add(new WeightedValue(4, 1));
		data.add(new WeightedValue(4, 2));
		data.add(new WeightedValue(4, 1));
		data.add(new WeightedValue(5, 0.5));
		data.add(new WeightedValue(5, 0.5));
		data.add(new WeightedValue(5, 0.5));
		data.add(new WeightedValue(5, 0.5));
		return data;
	}
}
