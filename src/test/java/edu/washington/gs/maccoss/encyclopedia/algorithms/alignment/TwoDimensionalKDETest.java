package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter3d;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.MedianInterpolatorTest;
import gnu.trove.list.array.TFloatArrayList;

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
		TwoDimensionalKDE filter=new TwoDimensionalKDE(rts);
		
		Charter3d.plot(filter, filter.getXRange(), filter.getYRange(), filter.getResolution()/5);
		
		plot(rts, filter.trace());
	}

	
	public static void plot(ArrayList<XYPoint> rts, Function rtWarper) {
		TFloatArrayList deltas=new TFloatArrayList();
		ArrayList<XYPoint> removedRTs=new ArrayList<XYPoint>();
		ArrayList<XYPoint> selectedRTs=new ArrayList<XYPoint>();
		for (int i=0; i<rts.size(); i++) {
			XYPoint xyPoint=rts.get(i);
			float actualRT=(float)xyPoint.y;
			float modelRT=rtWarper.getYValue((float)xyPoint.x);
			float delta=actualRT-modelRT;
			deltas.add(delta);

			removedRTs.add(xyPoint);
		}
		
		XYTrace median2=new XYTrace(rtWarper.getKnots(), GraphType.line, "Retention Time Fit");
		XYTrace selectedTrace=new XYTrace(selectedRTs, GraphType.tinypoint, "Data Used In Fit");
		XYTrace trace=new XYTrace(removedRTs, GraphType.tinypoint, "Data Removed From Fit");
		
		Charter.launchChart("Library RT", "Actual RT", true, median2, selectedTrace, trace);
	}
}
