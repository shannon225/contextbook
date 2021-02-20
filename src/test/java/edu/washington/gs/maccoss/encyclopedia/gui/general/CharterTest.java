package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.Color;
import java.awt.Dimension;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class CharterTest {
	public static void main(String[] args) {
		File[] fs=new File("/Volumes/searle_ssd/malaria/UW_SCX_yeast/scx_mzmls").listFiles(new SimpleFilenameFilter(".txt"));
		Arrays.sort(fs);
		for (File f : fs) {
			System.out.println(f.getName());
			HashMap<Comparable, TFloatArrayList> hash=new HashMap<>();
			TableParser.parseSSV(f, new TableParserMuscle() {
				@Override
				public void processRow(Map<String, String> row) {
					Integer time=Math.round((Integer.parseInt(row.get("time"))/10.0f))*10;
					String s = row.get("count");
					if (s==null||s.length()==0) return;
					float point=Float.parseFloat(s);
					if (point>20.0f) point=20.0f;
					TFloatArrayList list=hash.get(time);
					if (list==null) {
						list=new TFloatArrayList();
						hash.put(time, list);
					}
					list.add(point);
				}
				@Override
				public void cleanup() {
				}
			});
			
			System.out.println("Parsed "+hash.size()+" values...");
			final ExtendedChartPanel times=Charter.getBoxplotChart(null, "Retention Time (min)", "Number of MS2s", hash);
			Charter.launchChart(times, "Number of MS2s", new Dimension(500, 200));
		}

	}
	
	
	public static void main2(String[] args) {
		File f=new File("/Users/searleb/Documents/chromatogram_library_manuscript/quant_replicates/cvs_by_mean.txt");
		
		final TFloatFloatHashMap hash=new TFloatFloatHashMap();
		TableParser.parseTSV(f, new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				hash.put(Log.log10(Float.parseFloat(row.get("mean"))), Float.parseFloat(row.get("cv")));
			}
			@Override
			public void cleanup() {
			}
		});
		
		System.out.println("Parsed "+hash.size()+" values...");
		XYTraceInterface dataset=new XYTrace(hash, GraphType.tinypoint, "Coefficient of Variance", new Color(0.0f, 0.0f, 1.0f, 0.1f), 0.1f);
		Charter.launchChart("Log10 Average Intensity", "Coefficient of Variance", false, dataset);
	}

	/**
	 * Methionine: http://mona.fiehnlab.ucdavis.edu/spectra/display/CE000452
	 * @param args
	 */
	public static void main3(String[] args) {
		double[] masses = new double[] { 104.052841, 133.031921, 150.038239, 150.046448, 150.04895, 150.050446,
				150.051987, 150.054855, 150.058304, 150.063202, 150.065689, 150.067154, 150.069916, 150.074142,
				299.242279 };
		float[] intensities = new float[] { 9.029358f, 41.774054f, 0.503758f, 0.418014f, 0.40345f, 0.539593f, 0.682148f,
				1.009503f, 100f, 0.635989f, 0.550531f, 0.604124f, 0.460001f, 0.433143f, 0.651692f};
		
		Charter.launchChart(new LibraryEntry("Massbank", new HashSet<>(), 150.0583, (byte)1, "M", 1, 3791.84f, 0, masses, intensities, new AminoAcidConstants()));
	}
}
