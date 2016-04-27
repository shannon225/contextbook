package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ProphetMixtureModel;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RunningMedianWarper;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Gaussian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.UnitDistribution;
import gnu.trove.list.array.TFloatArrayList;

public class RetentionTimeFilter {
	private final RunningMedianWarper rtWarper;
	private final ProphetMixtureModel model;
	
	public RetentionTimeFilter(ArrayList<XYPoint> rts) {
		Collections.sort(rts);
		int order=Math.max(3, Math.round(rts.size()/50f));
		
		// 5 iterations is sufficient
		RunningMedianWarper warper=new RunningMedianWarper(rts, order, true);
		Pair<RunningMedianWarper, ProphetMixtureModel> pair=runIteration(rts, warper);
		pair=runIteration(rts, pair.x);
		pair=runIteration(rts, pair.x);
		pair=runIteration(rts, pair.x);
		rtWarper=pair.x;
		model=pair.y;
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
			deltas.add(delta);
			
			float prob=getProbabilityFitsModel((float)xyPoint.y, (float)xyPoint.x);
			if (prob>=0.5f) {
				selectedRTs.add(xyPoint);
			} else {
				removedRTs.add(xyPoint);
			}
		}
		float[] deltaArray=deltas.toArray();
		ArrayList<XYPoint> histogram=PivotTableGenerator.createPivotTable(deltaArray);
		XYTrace histTrace=new XYTrace(histogram, GraphType.line, "Delta RT");
		
		ArrayList<XYPoint> positivePoints=new ArrayList<XYPoint>();
		ArrayList<XYPoint> negativePoints=new ArrayList<XYPoint>();
		for (XYPoint xyPoint : histogram) {
			double x=xyPoint.getX();
			positivePoints.add(new XYPoint(x, model.getPositive().getProbability(x)));
			negativePoints.add(new XYPoint(x, model.getNegative().getProbability(x)));
		}

		XYTrace posTrace=new XYTrace(positivePoints, GraphType.line, "Positive");
		XYTrace negTrace=new XYTrace(negativePoints, GraphType.line, "Negative");
		
		XYTrace median2=new XYTrace(rtWarper.getKnots(), GraphType.line, "Retention Time Fit");
		XYTrace selectedTrace=new XYTrace(selectedRTs, GraphType.tinypoint, "Data Used In Fit");
		XYTrace trace=new XYTrace(removedRTs, GraphType.tinypoint, "Data Removed From Fit");
		
		if (saveFileSeed.isPresent()) {
			String saveFilePrefix=saveFileSeed.get().getAbsolutePath();
			Charter.writeAsPDF(new File(saveFilePrefix+".delta_rt.pdf"), "Delta RT", "Count", true, negTrace, posTrace, histTrace);
			Charter.writeAsPDF(new File(saveFilePrefix+".rt_fit.pdf"), "Calculated RT", "Actual RT", true, trace, selectedTrace, median2);
		} else {
			Charter.launchChart("Delta RT", "Count", true, histTrace, posTrace, negTrace);
			Charter.launchChart("Calculated RT", "Actual RT", true, median2, selectedTrace, trace);
		}
	}
	
	public float getYValue(float xrt) {
		return rtWarper.getYValue(xrt);
	}
	
	public float getProbabilityFitsModel(float actualRT, float modelRT) {
		float delta=actualRT-getYValue(modelRT);
		float probability=model.getProbability(delta);
		return probability;
	}

	public Pair<RunningMedianWarper,ProphetMixtureModel> runIteration(ArrayList<XYPoint> rts, RunningMedianWarper warper) {
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
		Distribution positive=new Gaussian(median, iqr/1.35f, 0.5f);
		Distribution negative=new UnitDistribution(median, quarterMaxRange, 0.5f, min, max);
		
		ProphetMixtureModel model=new ProphetMixtureModel(positive, negative, true);
		model.train(deltaArray, 10);
		positive=model.getPositive();
		negative=model.getNegative();
		
		/*
		ArrayList<XYPoint> histogram=PivotTableGenerator.createPivotTable(deltaArray);
		XYTrace histTrace=new XYTrace(histogram, GraphType.line, "Delta RT");
		
		ArrayList<XYPoint> positivePoints=new ArrayList<XYPoint>();
		ArrayList<XYPoint> negativePoints=new ArrayList<XYPoint>();
		for (XYPoint xyPoint : histogram) {
			double x=xyPoint.getX();
			positivePoints.add(new XYPoint(x, positive.getProbability(x)));
			negativePoints.add(new XYPoint(x, negative.getProbability(x)));
		}

		XYTrace posTrace=new XYTrace(positivePoints, GraphType.line, "Positive");
		XYTrace negTrace=new XYTrace(negativePoints, GraphType.line, "Negative");
		
		Charter.launchChart("Delta RT", "Count", true, histTrace, posTrace, negTrace);
		*/
		
		ArrayList<XYPoint> selectedRTs=new ArrayList<XYPoint>();
		for (int i=0; i<rts.size(); i++) {
			XYPoint xyPoint=rts.get(i);
			float delta=(float)xyPoint.y-warper.getYValue((float)xyPoint.x);
			if (model.getProbability(delta)>=0.5f) {
				selectedRTs.add(xyPoint);
			}
		}

		int order=Math.max(3, Math.round(selectedRTs.size()/50f));
		warper=new RunningMedianWarper(selectedRTs, order, true);

		/*
		XYTrace median2=new XYTrace(warper.getKnots(), GraphType.line, "corrected");
		XYTrace selectedTrace=new XYTrace(selectedRTs, GraphType.tinypoint, "Selected");

		XYTrace trace=new XYTrace(rts, GraphType.tinypoint, "RT");
		Charter.launchChart("Calculated RT", "Actual RT", true, median2, selectedTrace, trace);
		*/
		return new Pair<RunningMedianWarper, ProphetMixtureModel>(warper, model);
	}
}
