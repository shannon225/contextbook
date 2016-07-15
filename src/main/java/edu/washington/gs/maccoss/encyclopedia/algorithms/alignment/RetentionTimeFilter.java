package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ProphetMixtureModel;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Gaussian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.UnitDistribution;
import gnu.trove.list.array.TFloatArrayList;

public class RetentionTimeFilter {
	public static final float maxDeltaForHistogram=10.0f; // in minutes
	public static final float rejectionPValue=0.05f;
	private final Function rtWarper;
	private final ProphetMixtureModel model;
	private final String xAxis,yAxis;
	
	public RetentionTimeFilter(ArrayList<XYPoint> rts) {
		this(rts, "Library RT", "Actual RT");
	}
	public RetentionTimeFilter(ArrayList<XYPoint> rts, String xAxis, String yAxis) {
		TwoDimensionalKDE twoDimKDE=new TwoDimensionalKDE(rts);
		rtWarper=twoDimKDE.trace();
		model=generateMixtureModel(rts, rtWarper);
		this.xAxis=xAxis;
		this.yAxis=yAxis;
	}
	
	public void plot(ArrayList<XYPoint> rts) {
		File seed=null;
		plot(rts, Optional.ofNullable(seed));
	}
	
	public void plot(ArrayList<XYPoint> rts, Optional<File> saveFileSeed) {
		TFloatArrayList deltas=new TFloatArrayList();
		ArrayList<XYPoint> removedRTs=new ArrayList<XYPoint>();
		ArrayList<XYPoint> selectedRTs=new ArrayList<XYPoint>();
		for (int i=0; i<rts.size(); i++) {
			XYPoint xyPoint=rts.get(i);
			float actualRT=(float)xyPoint.y;
			float modelRT=rtWarper.getYValue((float)xyPoint.x);
			float delta=actualRT-modelRT;
			
			if (delta>-maxDeltaForHistogram&&delta<maxDeltaForHistogram) {
				deltas.add(delta);
			}
			
			float prob=getProbabilityFitsModel((float)xyPoint.y, (float)xyPoint.x);
			if (prob>=rejectionPValue) {
				selectedRTs.add(xyPoint);
			} else {
				removedRTs.add(xyPoint);
			}
		}
		float[] deltaArray=deltas.toArray();
		Arrays.sort(deltaArray);
		int min=0; //Math.round(deltaArray.length*0.05f);
		int max=deltaArray.length-1; //Math.round(deltaArray.length*0.95f);
		//float[] truncatedDeltaArray=new float[max-min];
		//System.arraycopy(deltaArray, min, truncatedDeltaArray, 0, max-min);
		
		ArrayList<XYPoint> histogram=PivotTableGenerator.createPivotTable(deltaArray);
		XYTrace histTrace=new XYTrace(histogram, GraphType.area, "Delta RT", new Color(26, 198, 49), 3.0f);
		
		ArrayList<XYPoint> positivePoints=new ArrayList<XYPoint>();
		ArrayList<XYPoint> negativePoints=new ArrayList<XYPoint>();
		
		int numPoints=500;
		double range=deltaArray[max]-deltaArray[min];
		for (int i=0; i<numPoints; i++) {
			double x=deltaArray[min]+i*range/numPoints;
			positivePoints.add(new XYPoint(x, model.getPositive().getProbability(x)));
			negativePoints.add(new XYPoint(x, model.getNegative().getProbability(x)));
		}

		double histSum=0.0;
		for (XYPoint xyPoint : histogram) {
			histSum+=xyPoint.getY();
		}
		
		double distSum=0;
		for (XYPoint xyPoint : positivePoints) {
			distSum+=xyPoint.getY();
		}
		for (XYPoint xyPoint : negativePoints) {
			distSum+=xyPoint.getY();
		}
		double normalizer=distSum>0?(histSum*numPoints)/(distSum*histogram.size()):1.0;
		
		ArrayList<XYPoint> normPositivePoints=new ArrayList<XYPoint>();
		for (XYPoint xyPoint : positivePoints) {
			normPositivePoints.add(new XYPoint(xyPoint.x, xyPoint.y*normalizer));
		}
		positivePoints=normPositivePoints;
		
		ArrayList<XYPoint> normNegativePoints=new ArrayList<XYPoint>();
		for (XYPoint xyPoint : negativePoints) {
			normNegativePoints.add(new XYPoint(xyPoint.x, xyPoint.y*normalizer));
		}
		negativePoints=normNegativePoints;

		XYTrace posTrace=new XYTrace(positivePoints, GraphType.line, "Positive", Color.BLUE, 2.0f);
		XYTrace negTrace=new XYTrace(negativePoints, GraphType.line, "Negative", Color.RED, 2.0f);
		
		XYTrace median2=new XYTrace(rtWarper.getKnots(), GraphType.line, "Retention Time Fit", new Color(26, 198, 49, 100), 3.0f);
		XYTrace selectedTrace=new XYTrace(selectedRTs, GraphType.tinypoint, "Data Used In Fit", Color.BLUE, 3.0f);
		XYTrace trace=new XYTrace(removedRTs, GraphType.tinypoint, "Data Removed From Fit", Color.RED, 3.0f);
		
		if (saveFileSeed.isPresent()) {
			String saveFilePrefix=saveFileSeed.get().getAbsolutePath();
			Charter.writeAsPDF(new File(saveFilePrefix+".delta_rt.pdf"), "Delta RT", "Count", true, negTrace, posTrace, histTrace);
			Charter.writeAsPDF(new File(saveFilePrefix+".rt_fit.pdf"), xAxis, yAxis, true, median2, selectedTrace, trace);
		} else {
			Charter.launchChart("Delta RT", "Count", true, negTrace, posTrace, histTrace);
			Charter.launchChart(xAxis, yAxis, true, median2, selectedTrace, trace);
		}
	}
	
	public float getYValue(float xrt) {
		return rtWarper.getYValue(xrt);
	}
	public float getXValue(float yrt) {
		return rtWarper.getXValue(yrt);
	}
	
	public float getProbabilityFitsModel(float actualRT, float modelRT) {
		float delta=actualRT-getYValue(modelRT);
		float probability=model.getProbability(delta);
		return probability;
	}

	public ProphetMixtureModel generateMixtureModel(ArrayList<XYPoint> rts, Function warper) {
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
		//Arrays.sort(deltaArray);
		//min=deltaArray[Math.round(deltaArray.length*0.05f)];
		//max=deltaArray[Math.round(deltaArray.length*0.95f)];
		
		float median=QuickMedian.select(deltaArray, 0.5f);
		float iqr=QuickMedian.iqr(deltaArray);
		float quarterMaxRange=(max-min)/4.0f;
		Distribution positive=new Gaussian(median, iqr/1.35f, 0.5f);
		Distribution negative=new UnitDistribution(median, quarterMaxRange, 0.5f, min, max);
		
		ProphetMixtureModel model=new ProphetMixtureModel(positive, negative, true);
		model.train(deltaArray, 10);
		positive=model.getPositive();
		negative=model.getNegative();
		return model;
	}
}
