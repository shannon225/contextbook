package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.util.Collections;
import java.util.List;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearRegression;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RTProbabilityModel;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SimpleMixtureModel;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Gaussian;
import gnu.trove.list.array.TFloatArrayList;

public class RetentionTimeTargetDecoyFilter extends AbstractRetentionTimeFilter {
	private static final float rejectionPValue=0.25f;
	
	public static RetentionTimeTargetDecoyFilter getFilter(List<XYPoint> targetMatches, List<XYPoint> decoys) {
		return getFilter(targetMatches, decoys, RT_STRING, "Retention Time (min)");
	}
	public static RetentionTimeTargetDecoyFilter getFilter(List<XYPoint> targetMatches, List<XYPoint> decoys, String xAxis, String yAxis) {
		return getFilter(targetMatches, decoys, xAxis, yAxis, TwoDimensionalKDE.DEFAULT_RESOLUTION);
	}
	public static RetentionTimeTargetDecoyFilter getFilter(List<XYPoint> targetMatches, List<XYPoint> decoys, String xAxis, String yAxis, int resolution) {
		Function rtWarper;
		Optional<RTProbabilityModel> model;
		if (targetMatches.size()>20) {
			Logger.logLine("Enough data points ("+targetMatches.size()+") to perform KDE alignment.");
			TwoDimensionalKDE twoDimKDE=new TwoDimensionalKDE(targetMatches, resolution);
			rtWarper=twoDimKDE.trace();
			model=Optional.of(generateMixtureModel(targetMatches, decoys, rtWarper));
		} else {
			if (targetMatches.size()<=1) {
				Logger.errorLine("Not enough data points ("+targetMatches.size()+") to perform KDE alignment, forced to use one-to-one mapping!");
				rtWarper=new LinearRegression(new float[] {0, 1}, new float[] {0, 1});
				model=Optional.empty();
			} else {
				Logger.errorLine("Not enough data points ("+targetMatches.size()+") to perform KDE alignment, forced to use linear regression!");
				rtWarper=new LinearRegression(targetMatches);
				model=Optional.of(generateMixtureModel(targetMatches, decoys, rtWarper));
			}
		}
		return new RetentionTimeTargetDecoyFilter(rtWarper, model, xAxis, yAxis);
	}
	
	private RetentionTimeTargetDecoyFilter(Function rtWarper, Optional<RTProbabilityModel> model, String xAxis, String yAxis) {
		super(rtWarper, model, xAxis, yAxis);
	}
	
	@Override
	public float getRejectionPValue() {
		return rejectionPValue;
	}

	public static RTProbabilityModel generateMixtureModel(List<XYPoint> targetMatches, List<XYPoint> decoys, Function warper) {
		TFloatArrayList decoyDeltas=new TFloatArrayList();
		for (int i=0; i<decoys.size(); i++) {
			XYPoint xyPoint=decoys.get(i);
			float delta=(float)xyPoint.y-warper.getYValue((float)xyPoint.x);
			decoyDeltas.add(delta);
		}

		float[] data=decoyDeltas.toArray();
		Gaussian negative=new Gaussian(General.mean(data), General.stdev(data), 0.5f);
		
		TFloatArrayList targetDeltas=new TFloatArrayList();
		for (int i=0; i<targetMatches.size(); i++) {
			XYPoint xyPoint=targetMatches.get(i);
			float delta=(float)xyPoint.y-warper.getYValue((float)xyPoint.x);
			targetDeltas.add(delta);
		}
		targetDeltas.sort();

		float[] targets=targetDeltas.toArray();
		Gaussian positive=new Gaussian(General.mean(targets), General.stdev(targets), 0.5f);
		Logger.logLine("Starting with mixture model from pos: "+positive.toString()+" and neg: "+negative.toString()+" from "+targets.length+" data points");
		
		// remove the worst 5% of the data
		int numToRemove=targetMatches.size()/20;
		Logger.logLine("Removing the worst "+numToRemove+" of "+targetMatches.size()+" retention time matches.");
		for (int i = 0; i < numToRemove; i++) {
			if (Math.abs(targetDeltas.get(0))>Math.abs(targetDeltas.get(targetDeltas.size()-1))) {
				targetDeltas.removeAt(0);
			} else {
				targetDeltas.removeAt(targetDeltas.size()-1);
			}
		}
		
		targets=targetDeltas.toArray();
		positive=new Gaussian(General.mean(targets), General.stdev(targets), 0.5f);
		Logger.logLine("Generating mixture model from pos: "+positive.toString()+" and neg: "+negative.toString()+" from "+targets.length+" data points");
		SimpleMixtureModel model=new SimpleMixtureModel(positive, negative);
		
		return model;
	}
}
