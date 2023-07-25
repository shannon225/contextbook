package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Map;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.MedianInterpolatorTest;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ProphetMixtureModel;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SimpleMixtureModel;
import gnu.trove.list.array.TFloatArrayList;
import junit.framework.TestCase;

public class RetentionTimeTargetDecoyFilterTest {
	//public void testErrors() {
	public static void main(String[] args) {
		
		InputStream is=MedianInterpolatorTest.class.getResourceAsStream("/deltart/deltart_example.csv");
		
		final ArrayList<XYPoint> targets=new ArrayList<XYPoint>();
		final ArrayList<XYPoint> decoys=new ArrayList<XYPoint>();
		
		TableParser.parseCSV(is, new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				float library=Float.parseFloat(row.get("library"));
				float actual=Float.parseFloat(row.get("actual"));
				boolean td=row.get("istarget").equalsIgnoreCase("true");
				if (td) {
					targets.add(new XYPoint(library, actual));
				} else {
					decoys.add(new XYPoint(library, actual));
				}
			}
			
			@Override
			public void cleanup() {
			}
		});
		
		System.out.println(targets.size()+","+decoys.size());
		
		RetentionTimeTargetDecoyFilter filter=RetentionTimeTargetDecoyFilter.getFilter(targets, decoys);
		SimpleMixtureModel model=(SimpleMixtureModel)filter.getModel().get();

		Charter.launchChart(Charter.getChart(new Range(-50, 50), model.getPositive(), model.getNegative()), "prophet");
		
		TFloatArrayList xs=new TFloatArrayList();
		TFloatArrayList ys=new TFloatArrayList();
		for (float x = -50f; x < 50f; x=x+0.1f) {
			float y=model.getProbability(0.0f, x);
			xs.add(x);
			ys.add(y);
		}
		Charter.launchChart(Charter.getChart("RT", "Prob", false, new XYTrace(xs.toArray(), ys.toArray(), GraphType.line, "Probs")), "prophet prob");
	}
}
