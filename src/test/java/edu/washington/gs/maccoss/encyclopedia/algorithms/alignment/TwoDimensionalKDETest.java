package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter3d;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.math.MedianInterpolatorTest;

public class TwoDimensionalKDETest {
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
		rts=new ArrayList<XYPoint>(rts.subList(0, 10000));
		System.out.println(rts.size());
		File f=new File("/Users/searleb/Downloads/clib/23aug2017_hela_serum_timecourse_pool_wide_001_170829031834.mzML.encyclopedia.txt.rt_fit.txt");
		//f=new File("/Users/searleb/Downloads/23aug2017_hela_serum_timecourse_pool_wide_001_170829031834.mzML.encyclopedia.txt.rt_fit.txt");
		rts=MedianInterpolatorTest.getData(f);
		
		RetentionTimeAlignmentInterface filter=new RetentionTimeFilter(rts);
		TwoDimensionalKDE kde=new TwoDimensionalKDE(rts);
		
		//Charter3d.plot(kde, kde.getXRange(), kde.getYRange(), kde.getResolution()/5);
		filter.plot(rts, Optional.ofNullable((File)f));
	}
}
