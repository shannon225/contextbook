package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearRegression;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ProphetMixtureModel;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RTProbabilityModel;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Gaussian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.UnitDistribution;
import gnu.trove.list.array.TFloatArrayList;

public class MassErrorFilter extends AbstractMassErrorFilter {

	public static MassErrorFilter getFilter(MassTolerance tolerance, int msLevel, ArrayList<XYPoint> rts) {
		Function rtWarper;
		Optional<RTProbabilityModel> model;
		if (rts.size()>20) {
			Logger.logLine("Enough data points ("+rts.size()+") to perform KDE alignment.");
			MassErrorKDE twoDimKDE=new MassErrorKDE(rts);
			rtWarper=twoDimKDE.trace();
			model=Optional.of(generateMixtureModel(rts, rtWarper));
		} else {
			if (rts.size()<=1) {
				Logger.errorLine("Not enough data points ("+rts.size()+") to perform KDE alignment, skipping correction!");
				rtWarper=new LinearRegression(new float[] {0, 1}, new float[] {0, 1});
				model=Optional.empty();
			} else {
				Logger.errorLine("Not enough data points ("+rts.size()+") to perform KDE alignment, forced to use linear regression!");
				rtWarper=new LinearRegression(rts);
				model=Optional.of(generateMixtureModel(rts, rtWarper));
			}
		}
		return new MassErrorFilter(tolerance, msLevel, rtWarper, model);
	}

	private MassErrorFilter(MassTolerance tolerance, int msLevel, Function rtWarper, Optional<RTProbabilityModel> model) {
		super(tolerance,msLevel, rtWarper, model);
	}

	public static ProphetMixtureModel generateMixtureModel(ArrayList<XYPoint> rts, Function warper) {
		TFloatArrayList deltas=new TFloatArrayList();
		float min=Float.MAX_VALUE;
		float max=-Float.MAX_VALUE;
		for (int i=0; i<rts.size(); i++) {
			XYPoint xyPoint=rts.get(i);
			float delta=getCorrectedMassError(warper, (float)xyPoint.x, (float)xyPoint.y);
			
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

	public static float getCorrectedMassError(Function rtWarper, float actualRT, float massError) {
		return massError-rtWarper.getYValue(actualRT);
	}
}
