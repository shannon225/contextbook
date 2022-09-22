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
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RTProbabilityModel;
import gnu.trove.list.array.TFloatArrayList;

public class AbstractMassErrorFilter implements MassErrorAlignmentInterface {

	protected static final String RT_STRING="Retention Time (min)";
	public static final float rejectionPValue=0.05f;
	protected final Function rtWarper;
	protected final Optional<RTProbabilityModel> model;
	protected final MassTolerance tolerance;
	protected final int msLevel;
	
	AbstractMassErrorFilter(MassTolerance tolerance, int msLevel, Function rtWarper, Optional<RTProbabilityModel> model) {
		this.tolerance=tolerance;
		this.msLevel=msLevel;
		this.rtWarper=rtWarper;
		this.model=model;
	}
	
	public Function getRtWarper() {
		return rtWarper;
	}
	@Override
	public List<MassErrorDataPoint> plot(ArrayList<XYPoint> rts, Optional<File> saveFileSeed) {
		String xAxis="Retention Time (Min)";
		String yAxis="Mass Error ("+tolerance.getUnits()+")";
		
		TFloatArrayList rtValues=new TFloatArrayList();
		TFloatArrayList removedMassErrors=new TFloatArrayList();
		TFloatArrayList selectedMassErrors=new TFloatArrayList();
		ArrayList<XYPoint> removedPeptides=new ArrayList<XYPoint>();
		ArrayList<XYPoint> selectedPeptides=new ArrayList<XYPoint>();
		for (int i=0; i<rts.size(); i++) {
			XYPoint xyPoint=rts.get(i);
			rtValues.add((float)xyPoint.x);
			
			float prob=getProbabilityFitsModel((float)xyPoint.x, (float)xyPoint.y);
			if (prob>=rejectionPValue) {
				selectedPeptides.add(xyPoint);
				selectedMassErrors.add((float)xyPoint.y);
			} else {
				removedPeptides.add(xyPoint);
				removedMassErrors.add((float)xyPoint.y);
			}
		}
		
		// build histogram
		ArrayList<XYPoint>[] histograms=PivotTableGenerator.createPivotTables(new float[][] {selectedMassErrors.toArray(), removedMassErrors.toArray()}, false);
		ArrayList<XYPoint> posHist=histograms[0];
		ArrayList<XYPoint> negHist=histograms[1];
	
		XYTraceInterface histTrace=new XYTrace(negHist, GraphType.area, "Removed", Color.red, 3.0f);
		XYTraceInterface posHistTrace=new XYTrace(posHist, GraphType.area, "Used In Fit", Color.blue, 3.0f);
		
		float alpha=Math.min(1.0f, 5000.0f/rts.size());
		XYTraceInterface median2=new XYTrace(rtWarper.getKnots(), GraphType.line, "Mass Error Fit", new Color(26, 198, 49, 100), 4.0f);
		XYTraceInterface selectedTrace=new XYTrace(selectedPeptides, GraphType.tinypoint, "Data Used In Fit", new Color(0f, 0f, 1f, alpha), 1.0f);
		XYTraceInterface trace=new XYTrace(removedPeptides, GraphType.tinypoint, "Data Removed From Fit", new Color(1f, 0f, 0f, alpha), 1.0f);
		
		if (saveFileSeed.isPresent()) {
			String saveFilePrefix=saveFileSeed.get().getAbsolutePath();
			Charter.writeAsPDF(new File(saveFilePrefix+".delta_ms"+msLevel+".pdf"), yAxis, "Number of Peptides", false, posHistTrace, histTrace);
			Charter.writeAsPDF(new File(saveFilePrefix+".ms"+msLevel+"_fit.pdf"), xAxis, yAxis, false, median2, selectedTrace, trace);
	
			try {
				final File file = new File(saveFilePrefix + ".ms"+msLevel+"_fit.txt");
				PrintWriter writer=new PrintWriter(file, "UTF-8");
				writer.println(xAxis+"\t"+yAxis+"\tpredictedMassError\tcorrectedAbsoluteMassError\tfitProb\tisDecoy\tsequence");
				
				for (int i=0; i<rts.size(); i++) {
					XYPoint xyPoint=rts.get(i);
					float predictedMassError=rtWarper.getYValue((float)xyPoint.x);
					float correctedAbsoluteMassError=getCorrectedMassError((float)xyPoint.x, (float)xyPoint.y);
	
					float prob=getProbabilityFitsModel((float)xyPoint.x, correctedAbsoluteMassError);
					
					if (xyPoint instanceof PeptideXYPoint) {
						PeptideXYPoint rtPoint=(PeptideXYPoint)xyPoint;
						writer.println(xyPoint.x+"\t"+xyPoint.y+"\t"+predictedMassError+"\t"+correctedAbsoluteMassError+"\t"+prob+"\t"+rtPoint.isDecoy()+"\t"+rtPoint.getPeptideModSeq());
					} else {
						writer.println(xyPoint.x+"\t"+xyPoint.y+"\t"+predictedMassError+"\t"+correctedAbsoluteMassError+"\t"+prob+"\t?\t?");
					}
				}
				writer.flush();
				writer.close();
	
				// Use a list that's backed by the file we just wrote; this allows
				// the points to fall out of memory until this list is accessed.
				return new LazyFileReadingMassErrorDataList(file, "UTF-8");
			} catch (IOException e) {
				Logger.errorLine("Error writing mass error mapping file.");
				Logger.errorException(e);
	
				return Collections.emptyList();
			}
		} else {
			Charter.launchChart(yAxis, "Number of Peptides", true, posHistTrace, histTrace);
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
	public float getProbabilityFitsModel(float actualRT, float massError) {
		float correctedMassError=getCorrectedMassError(actualRT, massError);
		
		return getProbability(actualRT, correctedMassError);
	}

	@Override
	public float getCorrectedMassError(float actualRT, float massError) {
		return Math.abs(massError-getYValue(actualRT));
	}

	public float getProbability(float actualRT, float correctedMassError) {
		if (model.isPresent()) {
			return model.get().getProbability(actualRT, correctedMassError);
		} else {
			return 1f;
		}
	}

	
	/**
	 * Allow fetching data that's been written to disk.
	 */
	class LazyFileReadingMassErrorDataList extends AbstractList<MassErrorDataPoint> {
		private final File file;
		private final String encoding;

		boolean isOpen = false;
		int size = 0;
		List<MassErrorDataPoint> data;

		public LazyFileReadingMassErrorDataList(File file, String encoding) {
			this.file = file;
			this.encoding = encoding;
		}

		@Override
		public MassErrorDataPoint get(int index) {
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
				final ImmutableList.Builder<MassErrorDataPoint> list = ImmutableList.builder();

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

						list.add(MassErrorDataPoint.of(x, y, pred, delta, prob, decoy, pepModSeq));
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