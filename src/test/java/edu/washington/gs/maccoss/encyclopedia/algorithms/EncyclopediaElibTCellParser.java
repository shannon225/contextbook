package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.awt.Dimension;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.zip.DataFormatException;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.stat.inference.TestUtils;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.CoefficientOfVariationCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.LibraryReportExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.SampleCoordinate;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroupInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.ExtendedChartPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.BenjaminiHochberg;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class EncyclopediaElibTCellParser {
	private static final float fdrThreshold = 0.01f;
	private static final float pvalueThreshold = 0.01f;
	public static String[] sampleNames=new String[] {"Naive", "Acute D3", "Chronic D3", "Chronic D7"};
	public static int[] tests=new int[] {3};
	public static int[] controls=new int[] {0, 1};
	public static HashMap<String, SampleCoordinate> sampleKey=new HashMap<>();

	public static void main(String[] args) throws IOException, SQLException, DataFormatException {
		loadMap();
		
		File file=new File("/Users/searleb/Documents/OSU/projects/yi/051622/dia/051622_Mouse_Tcell_DIA_quant_reports.elib");
		File proteinReportFile=new File("/Users/searleb/Documents/OSU/projects/yi/051622/dia/051622_Mouse_Tcell_DIA_quant_reports.elib.proteins.txt");
		File stub=new File(file.getParent(), "day7_specific_boxplots");
		FileUtils.deleteDirectory(stub);
		
		stub.mkdirs();
		File[] testDirs=new File[tests.length];
		for (int j = 0; j < tests.length; j++) {
			testDirs[j]=new File(stub, sampleNames[tests[j]]+"_boxplots");
			testDirs[j].mkdirs();

		}
		
		CoefficientOfVariationCalculator cvCalculator=new CoefficientOfVariationCalculator(sampleKey, sampleNames, 0.2f);

		LibraryFile library=new LibraryFile();
		library.openFile(file);
		
		ArrayList<String> accessions=new ArrayList<String>();
		TDoubleArrayList pvalues=new TDoubleArrayList();
		ArrayList<TFloatArrayList[]> datasets=new ArrayList<TFloatArrayList[]>();
		
		TableParser.parseTSV(proteinReportFile, new TableParserMuscle() {
			
			@Override
			public void processRow(Map<String, String> row) {
				String proteinAccession=row.get("Protein");
				TFloatArrayList[] rawdata=new TFloatArrayList[sampleNames.length];
				TDoubleArrayList[] data=new TDoubleArrayList[sampleNames.length];
				for (int i = 0; i < data.length; i++) {
					rawdata[i]=new TFloatArrayList();
					data[i]=new TDoubleArrayList();
				}
				for (Entry<String, SampleCoordinate> entry : sampleKey.entrySet()) {
					String value=row.get(entry.getKey());
					float rawIntensity=Float.parseFloat(value);
					float intensity=rawIntensity;
					if (intensity>1000) {
						intensity=Log.log10(intensity);
					} else {
						intensity=3;
					}
					data[entry.getValue().getSampleIndex()].add(intensity);
					rawdata[entry.getValue().getSampleIndex()].add(rawIntensity);
				}
				
				ArrayList<double[]> dataset=new ArrayList<>();
				
				for (int i = 0; i < data.length; i++) {
					dataset.add(data[i].toArray());
				}
				accessions.add(proteinAccession);
				pvalues.add(TestUtils.oneWayAnovaPValue(dataset));
				datasets.add(rawdata);
			}
			
			@Override
			public void cleanup() {
			}
		});
		
		double[] fdrs=BenjaminiHochberg.calculateAdjustedPValues(pvalues.toArray());
		for (int i = 0; i < fdrs.length; i++) {
			if (fdrs[i]>fdrThreshold) {
				continue;
			}
			
			TFloatArrayList[] rawdata=datasets.get(i);
			
			double[] worstPValues=new double[tests.length];
			byte[] direction=new byte[tests.length];
			for (int t=0; t<tests.length; t++) {
				int testIndex=tests[t];

				TFloatArrayList test=rawdata[testIndex];
				float avgTest=General.mean(test.toArray());
				
				for (int j=0; j<controls.length; j++) {
					int controlIndex=controls[j];
					double pvalue=TestUtils.tTest(log10(General.toDoubleArray(test.toArray())), log10(General.toDoubleArray(rawdata[controlIndex].toArray())));
					worstPValues[t]=Math.max(worstPValues[t], pvalue);
					float avgControl=General.mean(rawdata[controlIndex].toArray());
					byte localDirection;
					if (avgTest>avgControl) {
						localDirection=1; 
					} else {
						localDirection=-1;
					}
					
					if (direction[t]==0) {
						direction[t]=localDirection;
					} else if (localDirection!=direction[t]) {
						direction[t]=-2; // error state;
					}
				}
			}
			
			for (int j = 0; j < worstPValues.length; j++) {
				if (worstPValues[j]<pvalueThreshold) {
					String accession=accessions.get(i).split(";")[0];
					
					String adjective="Flat";
					if (direction[j]==-1) {
						adjective="Down";
					} else if (direction[j]==1) {
						adjective="Up";
					} else if (direction[j]==-2) {
						// mixed
						continue;
					}
					
					System.out.println(sampleNames[tests[j]]+"\t"+adjective+"\t"+accession+"\t"+fdrs[i]+"\t"+worstPValues[j]);	
					ExtendedChartPanel panel=Charter.getBoxplotChart(accession, "", "Intensity", sampleNames, rawdata);
					File f=new File(testDirs[j], accession+".pdf");
					Charter.writeAsPDF(panel.getChart(), f, new Dimension(250, 150));
					
				}
			}
		}
	}

	public static double[] log10(double[] v) {
		double[] r=new double[v.length];
		for (int i=0; i<r.length; i++) {
			if (v[i]>1000) {
				r[i]=Log.log10(v[i]);
			} else {
				r[i]=3;
			}
		}
		return r;
	}

	public static void loadMap() {
		sampleKey.put("051622_Mouse_Tcell_DIA_naive_rep01.dia", new SampleCoordinate(0,0));
		sampleKey.put("051622_Mouse_Tcell_DIA_naive_rep02.dia", new SampleCoordinate(1,0));
		sampleKey.put("051622_Mouse_Tcell_DIA_naive_rep03.dia", new SampleCoordinate(2,0));
		sampleKey.put("051622_Mouse_Tcell_DIA_naive_rep04.dia", new SampleCoordinate(3,0));
		sampleKey.put("051622_Mouse_Tcell_DIA_naive_rep05.dia", new SampleCoordinate(4,0));
		sampleKey.put("051622_Mouse_Tcell_DIA_acuteD3_rep01.dia", new SampleCoordinate(0,1));
		sampleKey.put("051622_Mouse_Tcell_DIA_acuteD3_rep02.dia", new SampleCoordinate(1,1));
		sampleKey.put("051622_Mouse_Tcell_DIA_acuteD3_rep03.dia", new SampleCoordinate(2,1));
		sampleKey.put("051622_Mouse_Tcell_DIA_acuteD3_rep04.dia", new SampleCoordinate(3,1));
		sampleKey.put("051622_Mouse_Tcell_DIA_acuteD3_rep05.dia", new SampleCoordinate(4,1));
		sampleKey.put("051622_Mouse_Tcell_DIA_chronicD3_rep01.dia", new SampleCoordinate(0,2));
		sampleKey.put("051622_Mouse_Tcell_DIA_chronicD3_rep02.dia", new SampleCoordinate(1,2));
		sampleKey.put("051622_Mouse_Tcell_DIA_chronicD3_rep03.dia", new SampleCoordinate(2,2));
		sampleKey.put("051622_Mouse_Tcell_DIA_chronicD3_rep04.dia", new SampleCoordinate(3,2));
		sampleKey.put("051622_Mouse_Tcell_DIA_chronicD3_rep05.dia", new SampleCoordinate(4,2));
		sampleKey.put("051622_Mouse_Tcell_DIA_chronicD7_rep01.dia", new SampleCoordinate(0,3));
		sampleKey.put("051622_Mouse_Tcell_DIA_chronicD7_rep02.dia", new SampleCoordinate(1,3));
		sampleKey.put("051622_Mouse_Tcell_DIA_chronicD7_rep03.dia", new SampleCoordinate(2,3));
		sampleKey.put("051622_Mouse_Tcell_DIA_chronicD7_rep04.dia", new SampleCoordinate(3,3));
		sampleKey.put("051622_Mouse_Tcell_DIA_chronicD7_rep05.dia", new SampleCoordinate(4,3));
	}
}
