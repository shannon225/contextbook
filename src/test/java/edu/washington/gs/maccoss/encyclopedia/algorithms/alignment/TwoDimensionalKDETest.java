package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;

import org.jzy3d.plot3d.builder.Mapper;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter3d;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
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
		
		//RetentionTimeFilter filter=RetentionTimeFilter.getFilter(rts);
		//filter.plot(rts, Optional.ofNullable((File)null));
		
		//RTFitMixtureModel model=new RTFitMixtureModel(rts, filter.getRtWarper());
		
		//filter=BrudererRetentionTimeFilter.getFilter(rts);
		//filter.plot(rts, Optional.ofNullable((File)null));
		//filter.plot(rts, Optional.ofNullable((File)f));
		
		TwoDimensionalKDE kde=new TwoDimensionalKDE(data, 1000);

		Mapper mapper=new Mapper() {
			@Override
			public double f(double arg0, double arg1) {
				return kde.f(arg0, arg1);
			}
		};
		Charter3d.plot(mapper, 
				new org.jzy3d.maths.Range(kde.getXRange().getStart(), kde.getXRange().getStop()), 
				new org.jzy3d.maths.Range(kde.getYRange().getStart(), kde.getYRange().getStop()), 
				kde.getResolution()/5);
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
}
