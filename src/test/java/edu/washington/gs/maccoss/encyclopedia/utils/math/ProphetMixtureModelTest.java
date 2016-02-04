package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Gaussian;
import gnu.trove.list.array.TFloatArrayList;
import junit.framework.TestCase;

public class ProphetMixtureModelTest extends TestCase {
	public void testModel() throws Exception {
		ArrayList<XYPoint> rts=MedianInterpolatorTest.getPhosphoData();
		Collections.sort(rts);
		int order=Math.max(3, Math.round(rts.size()/50f));
		RunningMedianWarper warper=new RunningMedianWarper(rts, order, true);
		TFloatArrayList deltas=new TFloatArrayList();
		float min=Float.MAX_VALUE;
		float max=-Float.MAX_VALUE;
		for (int i=0; i<rts.size(); i++) {
			XYPoint xyPoint=rts.get(i);
			float delta=(float)xyPoint.y-warper.getYValue((float)xyPoint.x);
			deltas.add(delta);
			if (delta>max) max=delta;
			if (delta<min) min=delta;
		}
		float[] deltaArray=deltas.toArray();
		
		float median=QuickMedian.select(deltaArray, 0.5f);
		float iqr=QuickMedian.iqr(deltaArray);
		float quarterMaxRange=(max-min)/4.0f;
		
		Distribution positive=new Gaussian(median, iqr, 0.95f);
		Distribution negative=new Gaussian(median, quarterMaxRange, 0.05f);
		ProphetMixtureModel model=new ProphetMixtureModel(positive, negative, false);
		model.train(deltaArray, 10);
		positive=model.getPositive();
		negative=model.getNegative();
		assertTrue(positive.getPrior()>negative.getPrior());
		assertTrue(Math.abs(positive.getMean())<5);
		assertTrue(Math.abs(negative.getMean())<10);
		assertTrue(positive.getStdev()<negative.getStdev());
		
		ArrayList<XYPoint> histogram=PivotTableGenerator.createPivotTable(deltaArray);
		//XYTrace histTrace=new XYTrace(histogram, GraphType.line, "Delta RT");
		
		ArrayList<XYPoint> positivePoints=new ArrayList<XYPoint>();
		ArrayList<XYPoint> negativePoints=new ArrayList<XYPoint>();
		for (XYPoint xyPoint : histogram) {
			double x=xyPoint.getX();
			positivePoints.add(new XYPoint(x, positive.getProbability(x)));
			negativePoints.add(new XYPoint(x, negative.getProbability(x)));
		}

		//XYTrace posTrace=new XYTrace(positivePoints, GraphType.line, "Positive");
		//XYTrace negTrace=new XYTrace(negativePoints, GraphType.line, "Negative");
		
		//Charter.launchChart("Delta RT", "Count", true, histTrace, posTrace, negTrace);
		
		//Thread.sleep(Long.MAX_VALUE);
	}
}
