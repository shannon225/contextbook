package edu.washington.gs.maccoss.encyclopedia.algorithms.curve;

import java.awt.Dimension;
import java.io.File;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import junit.framework.TestCase;

public class DilutionCurveFitterTest extends TestCase {
	public static void main(String[] args) {
		testIncrementDensity();
	}

	public static void testIncrementDensity() {
		float[] values = new float[] { 6024.163657f, 6024.163657f, 1953.753912f, 6040.202599f, 6050.184087f,
				5356.342656f, 6050.184087f, 4285.658108f, 6040.202599f, 6040.202599f, 6091.450214f, 2199.785037f,
				3449.606169f, 5342.50356f, 5342.50356f, 2054.357067f, 3359.588526f, 1455.312016f, 1455.312016f,
				1455.312016f, 2983.780843f, 5658.44874f, 2983.780843f, 6091.450214f, 5658.44874f, 2983.780843f,
				5658.44874f, 3826.262371f, 3826.262371f, 1006.464941f, 6083.085275f, 6083.085275f, 4285.658108f,
				3826.262371f, 1955.175302f, 4572.590254f, 4572.590254f, 3396.991712f, 1955.175302f, 3396.991712f,
				3735.28178f, 3735.28178f, 2335.025086f, 2732.224165f, 2732.224165f, 5063.184808f, 6158.256267f,
				5063.184808f, 6158.256267f, 360.7100301f, 2571.830146f, 5063.184808f, 3365.164021f, 3338.876672f,
				3338.876672f, 5356.342656f };

		float[] assayRT=new float[Math.round(110*60)]; // N+W minutes in second increments
		for (int i = 0; i < assayRT.length; i++) {
			assayRT[i]=i/60f;
		}
		float[] assayDensity=new float[assayRT.length];
		
		float windowInMin=40f;
		for (int i = 0; i < values.length; i++) {
			assayDensity=DilutionCurveFitter.incrementDensity(values[i], windowInMin, assayDensity);
		}
		
		XYTrace trace=new XYTrace(assayRT, assayDensity, GraphType.area, "Scheduling density");
		ChartPanel panel=Charter.getChart("Retention Time (min)", "Number of Peptides", true, trace);
		Charter.launchChart(panel, "Assay Density");
		//Charter.writeAsPDF(panel.getChart(), new File(outputDirectory, "assay_density.pdf"), new Dimension(600, 300));
	}
}
