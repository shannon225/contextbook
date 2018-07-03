package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.Color;
import java.io.File;
import java.util.Map;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class CharterTest {
	public static void main(String[] args) {
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

}
