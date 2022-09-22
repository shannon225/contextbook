package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.MedianInterpolatorTest;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.CosineGaussian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution;
import junit.framework.TestCase;

public class TwoDimensionalKDETest extends TestCase {
	public static void main(String[] args) {
		ArrayList<XYPoint> rts=MedianInterpolatorTest.getSyntheticData();
		rts.addAll(MedianInterpolatorTest.getLowNoiseData());
		rts.addAll(MedianInterpolatorTest.getLowNoiseData());
		rts.addAll(MedianInterpolatorTest.getLowNoiseData());
		rts.addAll(MedianInterpolatorTest.getLowNoiseData());
		rts.addAll(MedianInterpolatorTest.getLowNoiseData());
		//rts=MedianInterpolatorTest.getSyntheticData();
		rts=MedianInterpolatorTest.getPhosphoData();
		//rts=MedianInterpolatorTest.getCleanData();
		//rts=new ArrayList<XYPoint>(rts.subList(0, 10000));
		//System.out.println(rts.size());
		//File f=new File("/Users/searleb/Downloads/23aug2017_hela_serum_timecourse_pool_wide_001_170829031834.mzML.encyclopedia.txt.first.rt_fit.txt");
		//f=new File("/Users/searleb/Downloads/23aug2017_hela_serum_timecourse_pool_wide_001_170829031834.mzML.encyclopedia.txt.rt_fit.txt");
		//rts=MedianInterpolatorTest.getData(f);
		//File f=new File("/Users/searleb/Documents/phospho_localization/data/110515_bcs_hela_phospho_starved_20mz_500_900.dia.encyclopedia.txt.rt_fit.txt");
		//File f=new File("/Users/searleb/Documents/chromatogram_library_manuscript/hela_window_size/2018may14_hela_window_size_test_BCS_hela_wide_400_1200_1.mzML.encyclopedia.txt.first.rt_fit.txt");
		//rts=MedianInterpolatorTest.getData(f, 1f);
		
		ArrayList<XYPoint> data=new ArrayList();
		data=rts;
		
		RetentionTimeFilter filter=RetentionTimeFilter.getFilter(rts);
		filter.plot(rts, Optional.ofNullable((File)null));
		
		//RTFitMixtureModel model=new RTFitMixtureModel(rts, filter.getRtWarper());
		
		//filter=BrudererRetentionTimeFilter.getFilter(rts);
		//filter.plot(rts, Optional.ofNullable((File)null));
		//filter.plot(rts, Optional.ofNullable((File)f));
		
		TwoDimensionalKDE kde=new TwoDimensionalKDE(data, 1000);

//		Mapper mapper=new Mapper() {
//			@Override
//			public double f(double arg0, double arg1) {
//				return kde.f(arg0, arg1);
//			}
//		};
//		Charter3d.plot(mapper, 
//				new org.jzy3d.maths.Range(kde.getXRange().getStart(), kde.getXRange().getStop()), 
//				new org.jzy3d.maths.Range(kde.getYRange().getStart(), kde.getYRange().getStop()), 
//				kde.getResolution()/5);
	}
	
	public void testStamp() {
		Distribution dist=new CosineGaussian(0.0, 5, 10);
		float[][] stamp=TwoDimensionalKDE.getStamp(dist);
		for (int i = 0; i < stamp.length; i++) {
			for (int j = 0; j < stamp[i].length; j++) {
				assertTrue(stamp[i][j]>=0.0f);
				assertTrue(stamp[i][j]<1.0f);
			}
		}
	}

	public void testStackOverflow() {
		int resolution = 10000;

		ArrayList<XYPoint> points = new ArrayList<>();
		final int max = 900;
		for (int i = 0; i< max; i++) {
			points.add(new XYPoint(i, i));
		}

		final TwoDimensionalKDE twoDimensionalKDE = new TwoDimensionalKDE(points, resolution);
		try {
			twoDimensionalKDE.traceNorthEast(0, 0, new ArrayList<>());
		} catch (StackOverflowError soe) {
			fail("Reproduced stack overflow in TwoDimensionalKDE!");
		}
	}

	public void testHistogram() {
		ArrayList<XYPoint> points = new ArrayList<>();
		final int max = 900;
		for (int i=0; i < max; i++) {
			points.add(new XYPoint(i, i-(i%7)));
		}
		int resolution = 10;
		final TwoDimensionalKDE twoDimensionalKDE = new TwoDimensionalKDE(points, resolution);
		final float[][] histogram = twoDimensionalKDE.getTwoDimensionalHistogram();

		float[][] expected = new float[][] {{78.5398f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
				{9.424778f, 147.6548f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
				{0.0f, 6.2831855f, 150.7964f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
				{0.0f, 0.0f, 3.1415927f, 153.938f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
				{0.0f, 0.0f, 0.0f, 0.0f, 153.938f, 3.1415927f, 0.0f, 0.0f, 0.0f, 0.0f},
				{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 157.0796f, 0.0f, 0.0f, 0.0f, 0.0f},
				{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 4.712389f, 152.3672f, 0.0f, 0.0f, 0.0f},
				{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.5707964f, 153.938f, 1.5707964f, 0.0f},
				{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 152.3672f, 4.712389f},
				{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 78.5398f}};

		for (int i=0; i<expected.length; i++) {
			for (int j=0; j<expected[i].length; j++) {
				assertEquals(expected[i][j], histogram[i][j], 0.0000005f);
			}
		}
	}

	public void testTraceNorthEast() {
		ArrayList<XYPoint> points = new ArrayList<>();
		final int max = 900;
		for (int i=0; i < max; i++) {
			points.add(new XYPoint(i+i%13,  (i-max/2)*(i-max/2)-(i%7)));
		}
		int resolution = 11;
		final TwoDimensionalKDE twoDimensionalKDE = new TwoDimensionalKDE(points, resolution);
		final Function trace = twoDimensionalKDE.trace();

		// spot-check some values
		assertEquals(547.4994f, trace.getXValue(600), 0.000005f);
		assertEquals(548.3961f, trace.getXValue(800), 0.000005f);
		assertEquals(548.4365f, trace.getXValue(809), 0.000005f);

		assertEquals(12308.685f, trace.getYValue(600), 0.000005f);
		assertEquals(97413.055f, trace.getYValue(800), 0.000005f);
		assertEquals(99420.234f, trace.getYValue(809), 0.000005f);

	}

	public void testTraceSouthWest() {
		ArrayList<XYPoint> points = new ArrayList<>();
		final int max = 900;
		for (int i=0; i < max; i++) {
			points.add(new XYPoint(i*Math.sin(i), -i*Math.cos(i)));
		}
		int resolution = 11;
		final TwoDimensionalKDE twoDimensionalKDE = new TwoDimensionalKDE(points, resolution);
		final Function trace = twoDimensionalKDE.trace();

		// spot-check some values
		assertEquals(410.02866f, trace.getXValue(600), 0.000005f);
		assertEquals(786.4395f, trace.getXValue(800), 0.000005f);
		assertEquals(795.42413f, trace.getXValue(809), 0.000005f);

		assertEquals(709.7225f, trace.getYValue(600), 0.000005f);
		assertEquals(813.58367f, trace.getYValue(800), 0.000005f);
		assertEquals(822.5991f, trace.getYValue(809), 0.000005f);

	}
}
