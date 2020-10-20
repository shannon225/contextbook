package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.io.File;
import java.util.ArrayList;
import java.util.Map;

import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import junit.framework.TestCase;

public class LinearDiscriminantAnalysisTest {
	public static void main(String[] args) {
		File f=new File("/Users/searleb/Downloads/lda_test.txt");
		
		ArrayList<float[]> target=new ArrayList<>();
		ArrayList<float[]> decoy=new ArrayList<>();
		
		TableParser.parseTSV(f, new TableParserMuscle() {
			
			@Override
			public void processRow(Map<String, String> row) {
				int td=Integer.parseInt(row.get("TD"));
				
				float primary=Float.parseFloat(row.get("primary"));
				float xCorrModel=Float.parseFloat(row.get("xCorrModel"));
				float sumOfSquaredErrors=Float.parseFloat(row.get("sumOfSquaredErrors"));
				
				if (td>0) {
					target.add(new float[] {primary, xCorrModel, sumOfSquaredErrors});
				} else {
					decoy.add(new float[] {primary, xCorrModel, sumOfSquaredErrors});
				}
			}
			
			@Override
			public void cleanup() {
				// TODO Auto-generated method stub
				
			}
		});
		
		LinearDiscriminantAnalysis lda=LinearDiscriminantAnalysis.buildModel(target, decoy);
		double[] coefficients=lda.getCoefficients();
		for (int i = 0; i < coefficients.length; i++) {
			System.out.println(i+": "+coefficients[i]);
		}
		System.out.println("c: "+lda.getConstant());
		System.out.println();
		for (float[] data : target) {
			System.out.println("1\t"+lda.getScore(data));
		}
		for (float[] data : decoy) {
			System.out.println("-1\t"+lda.getScore(data));
		}
		
//		0: -0.013720960802214612
//		1: -0.1075904065402189
//		2: 0.9523161866448213
//		c: -0.2812703285341387
	}
}
