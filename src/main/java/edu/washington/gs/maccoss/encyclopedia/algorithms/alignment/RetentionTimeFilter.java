package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearRegression;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ProphetMixtureModel;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RTProbabilityModel;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Gaussian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.UnitDistribution;
import gnu.trove.list.array.TFloatArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

public class RetentionTimeFilter extends AbstractRetentionTimeFilter {
	private final int size;
	
	public static RetentionTimeFilter getFilter(List<XYPoint> rts) {
		return getFilter(rts, RT_STRING, "Retention Time (min)");
	}
	public static RetentionTimeFilter getFilter(List<XYPoint> rts, String xAxis, String yAxis) {
		return getFilter(rts, xAxis, yAxis, TwoDimensionalKDE.DEFAULT_RESOLUTION);
	}
	public static RetentionTimeFilter getFilter(List<XYPoint> rts, String xAxis, String yAxis, int resolution) {
		Function rtWarper;
		Optional<RTProbabilityModel> model;
		if (rts.size()>20) {
			Logger.logLine("Enough data points ("+rts.size()+") to perform KDE alignment.");
			TwoDimensionalKDE twoDimKDE=new TwoDimensionalKDE(rts, resolution);
			rtWarper=twoDimKDE.trace();
			model=Optional.of(generateMixtureModel(rts, rtWarper));
		} else {
			if (rts.size()<=1) {
				Logger.errorLine("Not enough data points ("+rts.size()+") to perform KDE alignment, forced to use one-to-one mapping!");
				rtWarper=new LinearRegression(new float[] {0, 1}, new float[] {0, 1});
				model=Optional.empty();
			} else {
				Logger.errorLine("Not enough data points ("+rts.size()+") to perform KDE alignment, forced to use linear regression!");
				rtWarper=new LinearRegression(rts);
				model=Optional.of(generateMixtureModel(rts, rtWarper));
			}
		}
		return new RetentionTimeFilter(rtWarper, model, xAxis, yAxis, rts.size());
	}
	
	private RetentionTimeFilter(Function rtWarper, Optional<RTProbabilityModel> model, String xAxis, String yAxis, int size) {
		super(rtWarper, model, xAxis, yAxis);
		this.size=size;
	}
	
	public int size() {
		return size;
	}

	public static ProphetMixtureModel generateMixtureModel(List<XYPoint> rts, Function warper) {
		TFloatArrayList deltas=new TFloatArrayList();
		float min=Float.MAX_VALUE;
		float max=-Float.MAX_VALUE;
		for (int i=0; i<rts.size(); i++) {
			XYPoint xyPoint=rts.get(i);
			float delta=getDelta(warper, (float)xyPoint.y, (float)xyPoint.x);
			deltas.add(delta);
			if (delta>max) max=delta;
			if (delta<min) min=delta;
		}
		float[] deltaArray=deltas.toArray();
		Arrays.sort(deltaArray);
		
		// assumes 1% FDR
		int lowerIndex=Math.min(Math.round(deltaArray.length*0.005f), deltaArray.length-1);
		min=deltaArray[lowerIndex];
		int upperIndex=Math.min(Math.round(deltaArray.length*0.995f), deltaArray.length-1);
		max=deltaArray[upperIndex];
		
		//Arrays.sort(deltaArray);
		//min=deltaArray[Math.round(deltaArray.length*0.05f)];
		//max=deltaArray[Math.round(deltaArray.length*0.95f)];
		
		//float median=QuickMedian.select(deltaArray, 0.5f);
		//float iqr=QuickMedian.iqr(deltaArray);
		float median=deltaArray[Math.min(Math.round(deltaArray.length*0.50f), deltaArray.length-1)];
		float iqr=deltaArray[Math.min(Math.round(deltaArray.length*0.75f), deltaArray.length-1)]-deltaArray[Math.min(Math.round(deltaArray.length*0.25f), deltaArray.length-1)];
		
		float quarterMaxRange=(max-min)/4.0f;
		Distribution positive=new Gaussian(median, iqr/1.35f, 0.5f);
		Distribution negative=new UnitDistribution(median, quarterMaxRange, 0.5f, min, max);
		
		ProphetMixtureModel model=new ProphetMixtureModel(positive, negative, true);
		model.train(deltaArray, 10);
		positive=model.getPositive();
		negative=model.getNegative();
		return model;
	}

	public static float getDelta(Function rtWarper, float actualRT, float modelRT) {
		float one=actualRT-rtWarper.getYValue(modelRT);
		float two=rtWarper.getXValue(actualRT)-modelRT;
		if (Math.abs(one)<Math.abs(two)) {
			return one;
		} else {
			return two;
		}
	}
}
