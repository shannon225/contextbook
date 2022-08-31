package edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.io.FileConcatenatorTest;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import junit.framework.TestCase;

public class RocPlotTest extends TestCase {
	public static void main(String[] args) {
		File f1=new File("/Users/searle.30/Downloads/scoring_tests/xcorr/23aug2017_hela_serum_timecourse_pool_dda_001.dia.scribe.features.txt");
		File f2=new File("/Users/searle.30/Downloads/scoring_tests/scribe/23aug2017_hela_serum_timecourse_pool_dda_001.dia.scribe.features.txt");
		File f3=new File("/Users/searle.30/Downloads/scoring_tests/tandem/23aug2017_hela_serum_timecourse_pool_dda_001.dia.scribe.features.txt");
		File f4=new File("/Users/searle.30/Downloads/scoring_tests/spectralangle/23aug2017_hela_serum_timecourse_pool_dda_001.dia.scribe.features.txt");
		
		RocPlot xCorrLib = getRocPlotFromFeatureFile(f1, "xCorrLib");
		RocPlot lnInvSumOfSquaredErrors = getRocPlotFromFeatureFile(f1, "lnInvSumOfSquaredErrors");
		RocPlot xTandem = getRocPlotFromFeatureFile(f1, "primary");
		RocPlot contrastAngle = getRocPlotFromFeatureFile(f1, "contrastAngle");
		RocPlot dotProduct = getRocPlotFromFeatureFile(f1, "DotProduct");
		
		XYTrace[] traces=new XYTrace[] {dotProduct.getTrace()}; //xCorrLib.getTrace(), lnInvSumOfSquaredErrors.getTrace(), contrastAngle.getTrace(), xTandem.getTrace()};
		Charter.launchChart("FDR", "Number of PSMs", true, traces);
	}
	
	
	public void testRocPlot() throws Exception {
		final File f=FileConcatenatorTest.writeTempFile("/featurefiles/truncated_single_feature.txt");
		final String scoreName="xCorrLib";
		
		RocPlot plot = getRocPlotFromFeatureFile(f, scoreName);
		
		XYTrace trace=plot.getTrace();
		Pair<double[], double[]> xy=trace.toArrays();
		double[] x=xy.x;
		double[] y=xy.y;

		for (int i = 1; i < x.length; i++) {
			assert(x[i]>x[i-1]);
		}
		for (int i = 1; i < y.length; i++) {
			assert(y[i]>=y[i-1]);
		}
	}


	private static RocPlot getRocPlotFromFeatureFile(final File f, final String scoreName) {
		final ArrayList<ScoredObject<Boolean>> scores=new ArrayList<>();
		TableParser.parseTSV(f, new TableParserMuscle() {
			
			@Override
			public void processRow(Map<String, String> row) {
				float score=Float.parseFloat(row.get(scoreName));
				boolean isTarget=Integer.parseInt(row.get("Label"))>0;
				scores.add(new ScoredObject<Boolean>(score, isTarget));
			}
			
			@Override
			public void cleanup() {
			}
		});
		
		Collections.sort(scores);
		Collections.reverse(scores);

		int numTargets=0;
		int numDecoys=0;
		RocPlot plot=new RocPlot(scoreName);
		for (ScoredObject<Boolean> score : scores) {
			if (score.y) {
				numTargets++;
			} else {
				numDecoys++;
			}

			float fdrValue=numTargets>0?numDecoys/(float)numTargets:1.0f;
			if (fdrValue>1.0f) fdrValue=1.0f;
			plot.addData(fdrValue, numTargets);
		}
		return plot;
	}

}
