package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.ImmutableList;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ProphetMixtureModel;
import gnu.trove.list.array.TFloatArrayList;

public class AbstractRetentionTimeFilter implements RetentionTimeAlignmentInterface {

	protected static final String RT_STRING="RT from Library";
	private static final String DELTA_RETENTION_TIME_STRING="Delta RT from Library (min)";
	public static final float maxDeltaForHistogram=10.0f;
	public static final float rejectionPValue=0.05f;
	protected final Function rtWarper;
	protected final Optional<ProphetMixtureModel> model;
	protected final String xAxis;
	protected final String yAxis;
	
	AbstractRetentionTimeFilter(Function rtWarper, Optional<ProphetMixtureModel> model, String xAxis, String yAxis) {
		this.rtWarper=rtWarper;
		this.model=model;
		this.xAxis=xAxis;
		this.yAxis=yAxis;
	}

	@Override
	public List<AlignmentDataPoint> plot(ArrayList<XYPoint> rts, Optional<File> saveFileSeed) {
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
		if (deltaArray.length==0) {
			Logger.errorLine("Sorry, not enough points to plot RT alignment");
			return Collections.emptyList();
		}
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
			if (model.isPresent()){
				positivePoints.add(new XYPoint(x, model.get().getPositive().getProbability(x)));
			}
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
		
		float alpha=Math.min(1.0f, 5000.0f/rts.size());
		XYTraceInterface median2=new XYTrace(rtWarper.getKnots(), GraphType.line, "Retention Time Fit", new Color(26, 198, 49, 100), 2.0f);
		XYTraceInterface selectedTrace=new XYTrace(selectedRTs, GraphType.tinypoint, "Data Used In Fit", new Color(0f, 0f, 1f, alpha), 1.0f);
		XYTraceInterface trace=new XYTrace(removedRTs, GraphType.tinypoint, "Data Removed From Fit", new Color(1f, 0f, 0f, alpha), 1.0f);

		//median2=new XYTrace(rtWarper.getKnots(), GraphType.line, "Retention Time Fit", new Color(255, 0, 0), 5.0f);
		//selectedTrace=new XYTrace(selectedRTs, GraphType.tinypoint, "Data Used In Fit", new Color(0f, 0f, 1f, alpha), 1.0f);
		//trace=new XYTrace(removedRTs, GraphType.tinypoint, "Data Removed From Fit", new Color(0f, 0f, 1f, alpha), 1.0f);
		
		if (saveFileSeed.isPresent()) {
			String saveFilePrefix=saveFileSeed.get().getAbsolutePath();
			Charter.writeAsPDF(new File(saveFilePrefix+".delta_rt.pdf"), DELTA_RETENTION_TIME_STRING, "Number of Peptides", false, posTrace, posHistTrace, histTrace);
			Charter.writeAsPDF(new File(saveFilePrefix+".rt_fit.pdf"), xAxis, yAxis, false, median2, selectedTrace, trace);
	
			try {
				final File file = new File(saveFilePrefix + ".rt_fit.txt");
				PrintWriter writer=new PrintWriter(file, "UTF-8");
				writer.println("library\tactual\twarpToActual\tdelta\tfitProb\tisDecoy\tsequence");
				
				for (int i=0; i<rts.size(); i++) {
					XYPoint xyPoint=rts.get(i);
					float modelRT=rtWarper.getYValue((float)xyPoint.x);
					float delta=getDelta((float)xyPoint.y, (float)xyPoint.x);
	
					float prob=getProbabilityFitsModel((float)xyPoint.y, (float)xyPoint.x);
					
					if (xyPoint instanceof RTRTPoint) {
						RTRTPoint rtPoint=(RTRTPoint)xyPoint;
						writer.println(xyPoint.x+"\t"+xyPoint.y+"\t"+modelRT+"\t"+delta+"\t"+prob+"\t"+rtPoint.isDecoy()+"\t"+rtPoint.getPeptideModSeq());
					} else {
						writer.println(xyPoint.x+"\t"+xyPoint.y+"\t"+modelRT+"\t"+delta+"\t"+prob+"\t?\t?");
					}
				}
				writer.flush();
				writer.close();
	
				// Use a list that's backed by the file we just wrote; this allows
				// the points to fall out of memory until this list is accessed.
				return new LazyFileReadingRtDataList(file, "UTF-8");
			} catch (IOException e) {
				Logger.errorLine("Error writing retention time mapping file.");
				Logger.errorException(e);
	
				return Collections.emptyList();
			}
		} else {
			Charter.launchChart("Delta RT", "Count", true, posTrace, posHistTrace, histTrace);
			Charter.launchChart(xAxis, yAxis, true, median2, selectedTrace, trace);
	
			return Collections.emptyList();
		}
	}

	@Override
	public float getYValue(float xrt) {
		return rtWarper.getYValue(xrt);
	}

	@Override
	public float getXValue(float yrt) {
		return rtWarper.getXValue(yrt);
	}

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

	@Override
	public float getProbabilityFitsModel(float delta) {
		if (model.isPresent()) {
			return model.get().getProbability(delta);
		} else {
			return 1f;
		}
	}

	
	/**
	 * Allow fetching data that's been written to disk.
	 */
	class LazyFileReadingRtDataList extends AbstractList<AlignmentDataPoint> {
		private final File file;
		private final String encoding;

		boolean isOpen = false;
		int size = 0;
		List<AlignmentDataPoint> data;

		public LazyFileReadingRtDataList(File file, String encoding) {
			this.file = file;
			this.encoding = encoding;
		}

		@Override
		public AlignmentDataPoint get(int index) {
			checkOpen();
			return data.get(index);
		}

		@Override
		public int size() {
			checkOpen();
			return size;
		}

		private synchronized void checkOpen() {
			if (!isOpen) {
				try {
					open();
				} catch (IOException e) {
					size = 0;
					data = Collections.emptyList();
				}
				isOpen = true;
			}
		}

		private synchronized void open() throws IOException {
			try (BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file), encoding))) {
				final ImmutableList.Builder<AlignmentDataPoint> list = ImmutableList.builder();

				String line;
				while (null != (line = reader.readLine())) {
					try {
						final String[] cells = line.split("\\t");
						if (cells.length != 7) {
							System.err.println("invalid line: " + line);
							continue;
						}

						float x = Float.parseFloat(cells[0]);
						float y = Float.parseFloat(cells[1]);
						float pred = Float.parseFloat(cells[2]);
						float delta = Float.parseFloat(cells[3]);
						float prob = Float.parseFloat(cells[4]);
						Boolean decoy = "?".equals(cells[5]) ? null : Boolean.parseBoolean(cells[5]);
						String pepModSeq = "?".equals(cells[6]) ? null : cells[6];

						list.add(AlignmentDataPoint.of(x, y, pred, delta, prob, decoy, pepModSeq));
					} catch (Exception e) {
						continue;
					}
				}

				data = list.build();
				size = data.size();
			}
		}
	}
}