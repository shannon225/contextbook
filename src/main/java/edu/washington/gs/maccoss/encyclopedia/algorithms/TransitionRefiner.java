package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import gnu.trove.list.array.TFloatArrayList;

public class TransitionRefiner {
	
	public static final float quantitativeCorrelationThreshold=0.9f;
	public static final float identificationCorrelationThreshold=0.75f;
	
	public static void main(String[] args) {
		ArrayList<float[]> chromatograms=new ArrayList<float[]>();
		float[] top=new float[] {78122.11f, 98378.266f, 142467.98f, 160690.98f, 222905.19f, 232847.81f, 307916.06f, 287985.97f, 354645.78f, 292500.47f, 363600.66f, 315389.8f, 347600.38f, 301699.44f, 294894.66f, 221306.75f, 186132.05f, 170257.3f, 146349.67f};
		float[] rts=new float[top.length];
		for (int i=0; i<rts.length; i++) {
			rts[i]=i;
		}
		chromatograms.add(top);
		chromatograms.add(new float[] {45688.473f, 52552.156f, 77305.87f, 75685.37f, 118698.836f, 108305.95f, 149139.72f, 148315.0f, 177143.22f, 185055.73f, 182493.88f, 157640.7f, 178510.28f, 137231.45f, 147941.67f, 128270.27f, 115709.18f, 144451.12f, 134215.73f});
		chromatograms.add(new float[] {12591.294f, 9210.12f, 12560.381f, 9796.7705f, 16568.979f, 13401.578f, 23627.754f, 31735.346f, 32540.912f, 44866.586f, 31987.303f, 26937.012f, 22891.26f, 14829.437f, 14344.313f, 11883.415f, 9714.611f, 11929.409f, 8152.291f});
		chromatograms.add(new float[] {369690.44f, 456645.72f, 651469.06f, 798680.9f, 1050497.0f, 1078652.0f, 1415159.4f, 1470529.2f, 1595412.0f, 1727540.2f, 1660950.8f, 1535321.1f, 1659728.6f, 1265447.8f, 1308786.6f, 952438.1f, 909168.4f, 707203.4f, 555946.06f});
		chromatograms.add(new float[] {12269.271f, 22352.746f, 25379.63f, 35125.26f, 40516.137f, 50424.46f, 58947.516f, 62769.434f, 62271.68f, 68998.4f, 79224.35f, 71739.13f, 84123.71f, 71447.32f, 78636.96f, 65355.36f, 63741.344f, 69759.18f, 71029.445f});
		chromatograms.add(new float[] {0.0f, 0.0f, 0.0f, 2272.649f, 2853.767f, 2879.7114f, 2626.8823f, 0.0f, 3799.889f, 7058.896f, 9528.616f, 0.0f, 2981.3376f, 1523.2601f, 6134.175f, 955.7171f, 0.0f, 0.0f, 2771.0017f});
		chromatograms.add(new float[] {3496.0784f, 3709.9358f, 9963.435f, 7364.1235f, 14094.167f, 10891.73f, 20391.303f, 6583.6714f, 20247.096f, 19192.396f, 17108.42f, 12707.958f, 21243.793f, 9334.771f, 15843.864f, 3872.4568f, 11016.68f, 4989.8164f, 7633.1284f});
		chromatograms.add(new float[] {88673.75f, 108412.02f, 134443.86f, 186275.72f, 229804.58f, 255614.5f, 297812.0f, 274456.28f, 308098.8f, 310620.28f, 303670.72f, 276719.75f, 319922.88f, 289453.9f, 259773.83f, 232285.22f, 194102.6f, 149186.31f, 134754.16f});
		chromatograms.add(new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 10589.013f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f});
		chromatograms.add(new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 11368.357f, 0.0f, 0.0f, 23255.867f, 0.0f, 0.0f, 0.0f});
		
		Triplet<float[], float[], Range> trio=identifyTransitions("AAPQS[+80.0]PSVPK", chromatograms, rts, true);
		float[] correlations=trio.x;
		float[] integrations=trio.y;
		for (int i=0; i<integrations.length; i++) {
			System.out.println(correlations[i]+"\t"+integrations[i]);
		}
	}

	public static Triplet<float[], float[], Range> identifyTransitions(String peptideModSeq, ArrayList<float[]> chromatograms, float[] retentionTimes) {
		return identifyTransitions(peptideModSeq, chromatograms, retentionTimes, false);
	}
	public static Triplet<float[], float[], Range> identifyTransitions(String peptideModSeq, ArrayList<float[]> chromatograms, float[] retentionTimes, boolean plot) {
		if (chromatograms.size()==0) return new Triplet<float[], float[], Range>(new float[0], new float[0], new Range(retentionTimes[0], retentionTimes[retentionTimes.length-1]));
		
		ArrayList<float[]> normalizedChromatograms=normalize(chromatograms);
		
		float[] medianChromatogram=new float[chromatograms.get(0).length];
		float maxMedian=0.0f;
		for (int i=0; i<medianChromatogram.length; i++) {
			TFloatArrayList list=new TFloatArrayList();
			for (float[] chromatogram : normalizedChromatograms) {
				if (chromatogram.length>i) { 
					list.add(chromatogram[i]);
				} else {
					list.add(0.0f);
				}
			}
			medianChromatogram[i]=QuickMedian.median(list.toArray());
			if (medianChromatogram[i]>maxMedian) {
				maxMedian=medianChromatogram[i];
			}
		}

		float threshold=maxMedian*0.01f; // 1% of max
		int firstData=-1;
		int lastData=-1;
		for (int i=0; i<medianChromatogram.length; i++) {
			if (medianChromatogram[i]>=threshold) {
				if (firstData==-1) {
					firstData=i;
				}
				if (i>lastData) {
					lastData=i;
				}
			}
		}
		
		int startIndex=firstData<=0?0:firstData-1;
		int stopIndex=lastData>=medianChromatogram.length-1?medianChromatogram.length-1:lastData+1;
		Range range=new Range(retentionTimes[startIndex], retentionTimes[stopIndex]);

		float medianMean=General.mean(medianChromatogram, startIndex, stopIndex);
		float[] correlationArray=new float[normalizedChromatograms.size()];
		float[] integrationArray=new float[correlationArray.length];
		for (int i=0; i<normalizedChromatograms.size(); i++) {
			float[] normalizedChromatogram=normalizedChromatograms.get(i);
			float fragmentMean=General.mean(normalizedChromatogram, startIndex, stopIndex);
			
			float medianDeltaSquareSum=0.0f;
			float fragmentDeltaSquareSum=0.0f;
			float deltaProductSum=0.0f;
			for (int j=startIndex; j<=stopIndex; j++) {
				float deltaMedian=medianChromatogram[j]-medianMean;
				float deltaFragment=normalizedChromatogram[j]-fragmentMean;
				medianDeltaSquareSum+=deltaMedian*deltaMedian;
				fragmentDeltaSquareSum+=deltaFragment*deltaFragment;
				deltaProductSum+=deltaMedian*deltaFragment;
			}
			// calculate correlation
			correlationArray[i]=deltaProductSum/((float)Math.sqrt(medianDeltaSquareSum*fragmentDeltaSquareSum));
			
			// calculate area
			float[] chromatogram=chromatograms.get(i);
			integrationArray[i]=0.0f;
			for (int j=startIndex+1; j<=stopIndex; j++) {
				float trapezoid=(retentionTimes[j]-retentionTimes[j-1])*(chromatogram[j-1]+chromatogram[j])/2.0f;
				integrationArray[i]+=trapezoid;
			}
		}
		
		if (plot) {
			XYTrace start=toBoundaries(firstData-1, "start");
			XYTrace stop=toBoundaries(lastData+1, "stop");

			HashMap<String, ChartPanel> panels=new HashMap<String, ChartPanel>();
			panels.put("unnormalized", getChart(chromatograms, correlationArray, start, stop));
			panels.put("normalized", getChart(normalizedChromatograms, correlationArray, start, stop));
			panels.put("median", Charter.getChart("scan", "intensity", false, toXYTrace(medianChromatogram, "median", null), start, stop));
			Charter.launchCharts(peptideModSeq+" chart", panels);
		}
		
		return new Triplet<float[], float[], Range>(correlationArray, integrationArray, range);
	}

	public static ChartPanel getChart(ArrayList<float[]> chromatograms, float[] correlationArray, XYTrace start, XYTrace stop) {
		XYTrace[] xytraces=new XYTrace[chromatograms.size()+2];
		for (int i=0; i<chromatograms.size(); i++) {
			float[] fs=chromatograms.get(i);
			
			Color c;
			if (correlationArray[i]>quantitativeCorrelationThreshold) {
				c=new Color(0, 205, 0);
			} else if (correlationArray[i]>identificationCorrelationThreshold) {
				c=new Color(255, 215, 0);
			} else {
				c=Color.red;
			}

			xytraces[i]=toXYTrace(fs, ""+i, c);
		}
		xytraces[xytraces.length-2]=start;
		xytraces[xytraces.length-1]=stop;
		ChartPanel panel=Charter.getChart("scan", "intensity", true, xytraces);
		return panel;
	}
	public static XYTrace toBoundaries(float f, String name) {
		ArrayList<XYPoint> points=new ArrayList<XYPoint>();
		points.add(new XYPoint(f, -0.1f));
		points.add(new XYPoint(f, 0.1f));
		return new XYTrace(points, GraphType.line, name);
	}

	public static XYTrace toXYTrace(float[] fs, String name, Color color) {
		ArrayList<XYPoint> points=new ArrayList<XYPoint>();
		for (int j=0; j<fs.length; j++) {
			XYPoint point=new XYPoint(j, fs[j]);
			points.add(point);
		}
		XYTrace trace=new XYTrace(points, GraphType.line, name, color);
		return trace;
	}

	public static ArrayList<float[]> normalize(ArrayList<float[]> chromatograms) {
		ArrayList<float[]> normalizedChromatograms=new ArrayList<float[]>();
		for (float[] fs : chromatograms) {
			normalizedChromatograms.add(General.normalize(fs));
		}
		return normalizedChromatograms;
	}
}
