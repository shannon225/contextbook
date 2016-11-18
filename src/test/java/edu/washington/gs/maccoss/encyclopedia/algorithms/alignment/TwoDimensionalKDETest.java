package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;

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
		File f=new File("/Users/searleb/Documents/school/projects/freezer/zeroDay/122715_bcs_hela_24mz_400_1000.dia.encyclopedia.txt.rt_fit.txt");
		//f=new File("/Users/searleb/Documents/school/projects/freezer/zeroDay/122715_bcs_hela_24mz_400_1000_dda.dia.encyclopedia.txt.rt_fit.txt");
		//rts=MedianInterpolatorTest.getData(f);
		
		RetentionTimeFilter filter=new RetentionTimeFilter(rts);
		TwoDimensionalKDE kde=new TwoDimensionalKDE(rts);
		
		Charter3d.plot(kde, kde.getXRange(), kde.getYRange(), kde.getResolution()/5);
		//filter.plot(rts, Optional.ofNullable((File)f));
	}
}
