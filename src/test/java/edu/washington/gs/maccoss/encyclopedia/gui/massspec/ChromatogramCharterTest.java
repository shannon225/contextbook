package edu.washington.gs.maccoss.encyclopedia.gui.massspec;

import java.awt.Color;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYDataset;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.ExtendedChartPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;

public class ChromatogramCharterTest { //extends TestCase {
	public static void main(String[] args) {
		Pair<Optional<ArrayList<XYTrace>>, Optional<ArrayList<XYTrace>>> pair=getFQAENLEHNKK();
		ChartPanel chartPanel=ChromatogramCharter.createChart(pair.x, pair.y);

		Charter.launchChart(chartPanel, "title");
	}
	
	public static Pair<Optional<ArrayList<XYTrace>>, Optional<ArrayList<XYTrace>>> getFQAENLEHNKK() {
		String peptide="FQAENLEHNKK";
		return new Pair<Optional<ArrayList<XYTrace>>, Optional<ArrayList<XYTrace>>>(getPeptidePrecursor(peptide), getPeptideFragment(peptide));
	}

	private static Optional<ArrayList<XYTrace>> getPeptidePrecursor(String peptide) {
		InputStream is=ChromatogramCharterTest.class.getResourceAsStream("/chromatograms/"+peptide+".precursor.txt");
		try {
			return Optional.of(getTraceData(is));
		} catch (Exception e) {
			Logger.logException(e);
			return Optional.empty();
		}
	}

	private static Optional<ArrayList<XYTrace>> getPeptideFragment(String peptide) {
		InputStream is=ChromatogramCharterTest.class.getResourceAsStream("/chromatograms/"+peptide+".fragment.txt");
		try {
			return Optional.of(getTraceData(is));
		} catch (Exception e) {
			Logger.logException(e);
			return Optional.empty();
		}
	}

	private static ArrayList<XYTrace> getTraceData(InputStream is) {
		final HashMap<String, ArrayList<XYPoint>> pointMap=new HashMap<>();
		TableParser.parseTSV(is, new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String rtTag="Retention Time (min)";
				String s = row.get(rtTag);
				if (s==null) {
					rtTag="Retention Time";
					s=row.get(rtTag);
				}
				double rt=Double.parseDouble(s);
				
				for (String column : row.keySet()) {
					if (!rtTag.equals(column)&&!"Row".equals(column)) {
						ArrayList<XYPoint> points=pointMap.get(column);
						if (points==null) {
							points=new ArrayList<>();
							pointMap.put(column, points);
						}
						points.add(new XYPoint(rt, Double.parseDouble(row.get(column))));
					}
				}
			}
			
			@Override
			public void cleanup() {
			}
		});
		
		ArrayList<XYTrace> traces=new ArrayList<>();
		for (Entry<String, ArrayList<XYPoint>> entry : pointMap.entrySet()) {
			String name=entry.getKey();
			Color col=RandomGenerator.randomColor(name.hashCode());
			traces.add(new XYTrace(entry.getValue(), GraphType.line, name, col, 2.0f));
		}
		return traces;
	}
}
