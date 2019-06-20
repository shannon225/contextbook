package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.Color;
import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.list.array.TDoubleArrayList;

public class MultiLineCharter {
	public static void main(String[] args) throws Exception {
		File f=new File("/Volumes/searle_ssd/malaria/malaria_dia_protein_ratios.txt");
		String delim="\t";
		
		BufferedReader in=new BufferedReader(new FileReader(f));

		int alpha=Math.round(10000f/2740);
		//int alpha=Math.round(10000f/1418);

		int count=0;
		String eachline;
		double[] x=null;
		ArrayList<XYTrace> traces=new ArrayList<>();
		
		while ((eachline=in.readLine())!=null) {
			count++;
			if (eachline.trim().length()==0) {
				continue;
			}
			
			if (count%100==0) {
				System.out.print(".");
			}
			
			boolean endedEarly=false;
			String[] data=eachline.split(delim, -1);
			TDoubleArrayList y=new TDoubleArrayList();
			for (int i=0; i<data.length; i++) {
				try {
					float value=Float.parseFloat(data[i]);
					if (Float.isFinite(value)&&value>0.0001) {
						y.add(Log.log10(value));
					} else {
						endedEarly=true;
						break;
					}
				} catch (Exception e) {
					// ignore
				}
			}
			
			if (x==null) {
				x=y.toArray();
				XYTrace trace=new XYTrace(x, x, GraphType.dashedline, "", new Color(0, 125, 255, 255), 3.0f);
				traces.add(trace);
			} else {
				double[] thisX=Arrays.copyOfRange(x, 0, y.size());
				XYTrace trace=new XYTrace(thisX, y.toArray(), GraphType.line, "", new Color(0, 0, 0, alpha), 1.0f);
				traces.add(trace);

				if (endedEarly) {
					XYTrace dot=new XYTrace(new double[] {thisX[thisX.length-1]}, new double[] {y.get(thisX.length-1)}, GraphType.point, "", new Color(0, 0, 0, alpha*5), 3.0f);
					traces.add(dot);
				}
			}
		}
		in.close();
		
		Charter.launchChart("X", "Y", false, new Dimension(300, 250), traces.toArray(new XYTrace[traces.size()]));
	}

}
