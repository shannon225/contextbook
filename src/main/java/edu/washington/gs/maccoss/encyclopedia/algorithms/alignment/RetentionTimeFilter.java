package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearRegression;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ProphetMixtureModel;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Gaussian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.UnitDistribution;
import gnu.trove.list.array.TFloatArrayList;

public class RetentionTimeFilter implements RetentionTimeAlignmentInterface {
	//private static final String RT_STRING="iRT from DDA Library";
	//private static final String DELTA_RETENTION_TIME_STRING="Delta RT from DDA Library (min)";
	//private static final String RT_STRING="RT from Chromatogram Library (min)";
	//private static final String DELTA_RETENTION_TIME_STRING="Delta RT from Chromatogram Library (min)";
	private static final String RT_STRING="RT from Library";
	private static final String DELTA_RETENTION_TIME_STRING="Delta RT from Library (min)";
	public static final float maxDeltaForHistogram=10.0f; // in minutes
	public static final float rejectionPValue=0.05f;
	private final Function rtWarper;
	private final ProphetMixtureModel model;
	private final String xAxis,yAxis;
	
	public RetentionTimeFilter(ArrayList<XYPoint> rts) {
		this(rts, RT_STRING, "Retention Time (min)");
	}
	
	public RetentionTimeFilter(ArrayList<XYPoint> rts, String xAxis, String yAxis) {
		if (rts.size()>20) {
			TwoDimensionalKDE twoDimKDE=new TwoDimensionalKDE(rts);
			rtWarper=twoDimKDE.trace();
		} else {
			if (rts.size()<=1) {
				rtWarper=new LinearRegression(new float[] {0, 1}, new float[] {0, 1});
			} else {
				rtWarper=new LinearRegression(rts);
			}
		}
		model=generateMixtureModel(rts, rtWarper);
		this.xAxis=xAxis;
		this.yAxis=yAxis;
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface#plot(java.util.ArrayList, java.util.Optional)
	 */
	@Override
	public void plot(ArrayList<XYPoint> rts, Optional<File> saveFileSeed) {
		TFloatArrayList deltas=new TFloatArrayList();
		ArrayList<XYPoint> removedRTs=new ArrayList<XYPoint>();
		ArrayList<XYPoint> selectedRTs=new ArrayList<XYPoint>();
		for (int i=0; i<rts.size(); i++) {
			XYPoint xyPoint=rts.get(i);
			float delta=getDelta((float)xyPoint.y, (float)xyPoint.x);
			
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
		ArrayList<XYPoint> posHist=new ArrayList<XYPoint>();
		ArrayList<XYPoint> negHist=new ArrayList<XYPoint>();
		
		for (XYPoint xyPoint : histogram) {
			float prob=getProbabilityFitsModel((float)xyPoint.x);
			if (prob>=rejectionPValue) {
				posHist.add(xyPoint);
				negHist.add(new XYPoint(xyPoint.x, 0.0));
			} else {
				negHist.add(xyPoint);
			}
		}

		XYTraceInterface histTrace=new XYTrace(negHist, GraphType.area, "Delta RT", Color.red, 3.0f);
		XYTraceInterface posHistTrace=new XYTrace(posHist, GraphType.area, "Delta RT", Color.blue, 3.0f);
		
		ArrayList<XYPoint> positivePoints=new ArrayList<XYPoint>();
		
		int numPoints=500;
		double range=deltaArray[max]-deltaArray[min];
		for (int i=0; i<numPoints; i++) {
			double x=deltaArray[min]+i*range/numPoints;
			positivePoints.add(new XYPoint(x, model.getPositive().getProbability(x)));
		}

		double histSum=0.0;
		for (XYPoint xyPoint : histogram) {
			histSum+=xyPoint.getY();
		}
		
		double distSum=0;
		for (XYPoint xyPoint : positivePoints) {
			distSum+=xyPoint.getY();
		}
		double normalizer=distSum>0?(histSum*numPoints)/(distSum*histogram.size()):1.0;
		
		ArrayList<XYPoint> normPositivePoints=new ArrayList<XYPoint>();
		for (XYPoint xyPoint : positivePoints) {
			normPositivePoints.add(new XYPoint(xyPoint.x, xyPoint.y*normalizer));
		}
		positivePoints=normPositivePoints;
		
		XYTraceInterface posTrace=new XYTrace(positivePoints, GraphType.line, "Positive", new Color(26, 198, 49, 100), 2.0f);
		
		XYTraceInterface median2=new XYTrace(rtWarper.getKnots(), GraphType.line, "Retention Time Fit", new Color(26, 198, 49, 100), 2.0f);
		XYTraceInterface selectedTrace=new XYTrace(selectedRTs, GraphType.tinypoint, "Data Used In Fit", Color.BLUE, 1.0f);
		XYTraceInterface trace=new XYTrace(removedRTs, GraphType.tinypoint, "Data Removed From Fit", Color.RED, 1.0f);
		
		if (saveFileSeed.isPresent()) {
			String saveFilePrefix=saveFileSeed.get().getAbsolutePath();
			Charter.writeAsPDF(new File(saveFilePrefix+".delta_rt.pdf"), DELTA_RETENTION_TIME_STRING, "Number of Peptides", false, posTrace, posHistTrace, histTrace);
			Charter.writeAsPDF(new File(saveFilePrefix+".rt_fit.pdf"), xAxis, yAxis, false, median2, selectedTrace, trace);

			try {
				PrintWriter writer=new PrintWriter(new File(saveFilePrefix+".rt_fit.txt"), "UTF-8");
				writer.println("library\tactual\twarpToActual\tdelta\tfitProb");
				for (int i=0; i<rts.size(); i++) {
					XYPoint xyPoint=rts.get(i);
					float modelRT=rtWarper.getYValue((float)xyPoint.x);
					float delta=getDelta((float)xyPoint.y, (float)xyPoint.x);

					float prob=getProbabilityFitsModel((float)xyPoint.y, (float)xyPoint.x);
					writer.println(xyPoint.x+"\t"+xyPoint.y+"\t"+modelRT+"\t"+delta+"\t"+prob);
				}
				writer.flush();
				writer.close();
			} catch (IOException e) {
				Logger.errorLine("Error writing retention time mapping file.");
				Logger.errorException(e);
			}
		} else {
			Charter.launchChart("Delta RT", "Count", true, posTrace, posHistTrace, histTrace);
			Charter.launchChart(xAxis, yAxis, true, median2, selectedTrace, trace);
		}
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface#getYValue(float)
	 */
	@Override
	public float getYValue(float xrt) {
		return rtWarper.getYValue(xrt);
	}
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface#getXValue(float)
	 */
	@Override
	public float getXValue(float yrt) {
		return rtWarper.getXValue(yrt);
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface#getProbabilityFitsModel(float, float)
	 */
	@Override
	public float getProbabilityFitsModel(float actualRT, float modelRT) {
		float delta=getDelta(actualRT, modelRT);
		
		return getProbabilityFitsModel(delta);
	}

	private float getDelta(float actualRT, float modelRT) {
		float one=actualRT-getYValue(modelRT);
		float two=getXValue(actualRT)-modelRT;
		if (Math.abs(one)<Math.abs(two)) {
			return one;
		} else {
			return two;
		}
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface#getProbabilityFitsModel(float)
	 */
	@Override
	public float getProbabilityFitsModel(float delta) {
		float probability=model.getProbability(delta);
		return probability;
	}

	public ProphetMixtureModel generateMixtureModel(ArrayList<XYPoint> rts, Function warper) {
		TFloatArrayList deltas=new TFloatArrayList();
		float min=Float.MAX_VALUE;
		float max=-Float.MAX_VALUE;
		for (int i=0; i<rts.size(); i++) {
			XYPoint xyPoint=rts.get(i);
			float delta=getDelta((float)xyPoint.y, (float)xyPoint.x);
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
}
