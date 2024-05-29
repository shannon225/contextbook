package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.util.List;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearRegression;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RTProbabilityModel;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SimpleMixtureModel;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Gaussian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.KDE;
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
		//KDE negative=new KDE(decoyDeltas.toArray(), 0.5);
		float[] data=decoyDeltas.toArray();
		Gaussian negative=new Gaussian(General.mean(data), General.stdev(data), 0.5f);
		
		TFloatArrayList targetDeltas=new TFloatArrayList();
		for (int i=0; i<targetMatches.size(); i++) {
			XYPoint xyPoint=targetMatches.get(i);
			float delta=(float)xyPoint.y-warper.getYValue((float)xyPoint.x);
			targetDeltas.add(delta);
		}
		float[] targets = targetDeltas.toArray();
		float bottom=QuickMedian.select(targets, 0.005f);
		float top=QuickMedian.select(targets, 0.995f);
		//KDE positive=new KDE(targets, 0.5); // ok if they get sorted/rearranged
		Gaussian positive=new Gaussian(General.mean(targets), General.stdev(targets), 0.5f);
		
		SimpleMixtureModel model=new SimpleMixtureModel(positive, negative, new Range(bottom, top));
		
		return model;
	}
}
