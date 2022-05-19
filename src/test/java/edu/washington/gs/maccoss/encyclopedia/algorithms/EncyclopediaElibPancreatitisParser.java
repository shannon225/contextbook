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

public class EncyclopediaElibPancreatitisParser {
	private static final float pvalueThreshold = 0.05f;
	public static String[] sampleNames=new String[] {"CP", "AP", "Cont 1", "Cont 2", "Frac"};
	public static int[] tests=new int[] {0, 1};
	public static int[] controls=new int[] {2, 3, 4};
	public static HashMap<String, SampleCoordinate> sampleKey=new HashMap<>();

	public static void main(String[] args) throws IOException, SQLException, DataFormatException {
		loadMap();
		
		File file=new File("/Users/searle.30/Documents/CCIC/maisam/032922_pancreatitis_grant_dataset/032922_pancreatitis_120_quant_reports.elib");
		File stub=new File(file.getParent(), "pancreatitis_120_boxplots");
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
		
		ArrayList<ProteinGroupInterface> proteinGroups=library.getProteinGroups();
		System.out.println(proteinGroups.size()+" total protein groups");
		File proteinReportFile=LibraryReportExtractor.extractMatrix(library, proteinGroups, true, Optional.of(cvCalculator));
		
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
			if (fdrs[i]>pvalueThreshold) {
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
		sampleKey.put("032922_pancreatitis_1ug_DIA_5156_AP_F_01.mzML", new SampleCoordinate(0,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5157_AP_M_01.mzML", new SampleCoordinate(1,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5159_AP_F_01.mzML", new SampleCoordinate(2,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5160_AP_F_01.mzML", new SampleCoordinate(3,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5164_AP_M_01.mzML", new SampleCoordinate(4,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5166_AP_F_01.mzML", new SampleCoordinate(5,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5167_AP_F_01.mzML", new SampleCoordinate(6,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5168_AP_F_01.mzML", new SampleCoordinate(7,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5170_AP_F_01.mzML", new SampleCoordinate(8,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5171_AP_M_01.mzML", new SampleCoordinate(9,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5172_AP_M_01.mzML", new SampleCoordinate(10,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5174_AP_F_01.mzML", new SampleCoordinate(11,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5177_AP_M_01.mzML", new SampleCoordinate(12,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5178_AP_F_01.mzML", new SampleCoordinate(13,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5180_AP_F_01.mzML", new SampleCoordinate(14,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5185_AP_M_01.mzML", new SampleCoordinate(15,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5186_AP_F_01.mzML", new SampleCoordinate(16,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5188_AP_F_01.mzML", new SampleCoordinate(17,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5205_AP_M_01.mzML", new SampleCoordinate(18,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5210_AP_F_01.mzML", new SampleCoordinate(19,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5212_AP_F_01.mzML", new SampleCoordinate(20,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5216_AP_F_01.mzML", new SampleCoordinate(21,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5221_AP_M_01.mzML", new SampleCoordinate(22,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5223_AP_M_01.mzML", new SampleCoordinate(23,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5231_AP_F_01.mzML", new SampleCoordinate(24,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5232_AP_F_01.mzML", new SampleCoordinate(25,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5234_AP_F_01.mzML", new SampleCoordinate(26,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_5235_AP_M_01.mzML", new SampleCoordinate(27,1));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC026_Cont_F_01.mzML", new SampleCoordinate(0,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC027_Cont_M_01.mzML", new SampleCoordinate(1,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC029_Cont_M_01.mzML", new SampleCoordinate(2,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_001_Cont_F_01.mzML", new SampleCoordinate(3,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_002_Cont_F_01.mzML", new SampleCoordinate(4,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_003_Cont_F_01.mzML", new SampleCoordinate(5,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_004_Cont_F_01.mzML", new SampleCoordinate(6,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_005_Cont_M_01.mzML", new SampleCoordinate(7,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_007_Cont_F_01.mzML", new SampleCoordinate(8,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_008_Cont_F_01.mzML", new SampleCoordinate(9,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_009_Cont_M_01.mzML", new SampleCoordinate(10,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_011_Cont_M_01.mzML", new SampleCoordinate(11,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_012_Cont_F_01.mzML", new SampleCoordinate(12,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_013_Cont_M_01.mzML", new SampleCoordinate(13,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_014_Cont_M_01.mzML", new SampleCoordinate(14,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_015_Cont_M_01.mzML", new SampleCoordinate(15,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_016_Cont_F_01.mzML", new SampleCoordinate(16,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_017_Cont_F_01.mzML", new SampleCoordinate(17,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_019_Cont_M_01.mzML", new SampleCoordinate(18,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_023_Cont_F_01.mzML", new SampleCoordinate(19,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_025_Cont_F_01.mzML", new SampleCoordinate(20,2));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC026_Cont2_F_01.mzML", new SampleCoordinate(0,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC027_Cont2_M_01.mzML", new SampleCoordinate(1,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC029_Cont2_M_01.mzML", new SampleCoordinate(2,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_001_Cont2_F_01.mzML", new SampleCoordinate(3,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_002_Cont2_F_01.mzML", new SampleCoordinate(4,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_004_Cont2_F_01.mzML", new SampleCoordinate(5,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_005_Cont2_M_01.mzML", new SampleCoordinate(6,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_007_Cont2_F_01.mzML", new SampleCoordinate(7,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_008_Cont2_F_01.mzML", new SampleCoordinate(8,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_009_Cont2_M_01.mzML", new SampleCoordinate(9,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_011_Cont2_M_01.mzML", new SampleCoordinate(10,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_012_Cont2_F_01.mzML", new SampleCoordinate(11,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_013_Cont2_M_01.mzML", new SampleCoordinate(12,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_014_Cont2_M_01.mzML", new SampleCoordinate(13,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_015_Cont2_M_01.mzML", new SampleCoordinate(14,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_016_Cont2_F_01.mzML", new SampleCoordinate(15,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_017_Cont2_F_01.mzML", new SampleCoordinate(16,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_019_Cont2_M_01.mzML", new SampleCoordinate(17,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_023_Cont2_F_01.mzML", new SampleCoordinate(18,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPC_025_Cont2_F_01.mzML", new SampleCoordinate(19,3));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_003_CP_M_01.mzML", new SampleCoordinate(0,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_004_CP_F_01.mzML", new SampleCoordinate(1,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_006_CP_F_01.mzML", new SampleCoordinate(2,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_009_CP_F_01.mzML", new SampleCoordinate(3,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_016_CP_M_01.mzML", new SampleCoordinate(4,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_022_CP_M_01.mzML", new SampleCoordinate(5,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_025_CP_F_01.mzML", new SampleCoordinate(6,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_027_CP_M_01.mzML", new SampleCoordinate(7,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_028_CP_F_01.mzML", new SampleCoordinate(8,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_031_CP_F_01.mzML", new SampleCoordinate(9,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_033_CP_M_01.mzML", new SampleCoordinate(10,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_034_CP_F_01.mzML", new SampleCoordinate(11,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_035_CP_F_01.mzML", new SampleCoordinate(12,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_037_CP_M_01.mzML", new SampleCoordinate(13,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_039_CP_F_01.mzML", new SampleCoordinate(14,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_041_CP_F_01.mzML", new SampleCoordinate(15,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_042_CP_F_01.mzML", new SampleCoordinate(16,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_045_CP_F_01.mzML", new SampleCoordinate(17,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_048_CP_M_01.mzML", new SampleCoordinate(18,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_051_CP_F_01.mzML", new SampleCoordinate(19,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_052_CP_F_01.mzML", new SampleCoordinate(20,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_053_CP_F_01.mzML", new SampleCoordinate(21,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_054_CP_M_01.mzML", new SampleCoordinate(22,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_057_CP_F_01.mzML", new SampleCoordinate(23,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_058_CP_M_01.mzML", new SampleCoordinate(24,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_061_CP_F_01.mzML", new SampleCoordinate(25,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_062_CP_M_01.mzML", new SampleCoordinate(26,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_063_CP_M_01.mzML", new SampleCoordinate(27,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_064_CP_M_01.mzML", new SampleCoordinate(28,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_065_CP_M_01.mzML", new SampleCoordinate(29,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_066_CP_M_01.mzML", new SampleCoordinate(30,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_068_CP_F_01.mzML", new SampleCoordinate(31,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_069_CP_M_01.mzML", new SampleCoordinate(32,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_070_CP_M_01.mzML", new SampleCoordinate(33,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_072_CP_F_01.mzML", new SampleCoordinate(34,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_073_CP_F_01.mzML", new SampleCoordinate(35,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_074_CP_M_01.mzML", new SampleCoordinate(36,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_075_CP_M_01.mzML", new SampleCoordinate(37,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_076_CP_F_01.mzML", new SampleCoordinate(38,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_078_CP_F_01.mzML", new SampleCoordinate(39,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_079_CP_M_01.mzML", new SampleCoordinate(40,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_080_CP_M_01.mzML", new SampleCoordinate(41,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_081_CP_M_01.mzML", new SampleCoordinate(42,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_083_CP_F_01.mzML", new SampleCoordinate(43,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_085_CP_F_01.mzML", new SampleCoordinate(44,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_088_CP_M_01.mzML", new SampleCoordinate(45,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_091_CP_F_01.mzML", new SampleCoordinate(46,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_092_CP_M_01.mzML", new SampleCoordinate(47,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_095_CP_F_01.mzML", new SampleCoordinate(48,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_2016_9510_096_CP_M_01.mzML", new SampleCoordinate(49,0));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF050_Frac_M_01.mzML", new SampleCoordinate(0,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF051_Frac_F_01.mzML", new SampleCoordinate(1,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF054_Frac_F_01.mzML", new SampleCoordinate(2,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF057_Frac_M_01.mzML", new SampleCoordinate(3,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF058_Frac_M_01.mzML", new SampleCoordinate(4,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF060_Frac_F_01.mzML", new SampleCoordinate(5,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_002_Frac_F_01.mzML", new SampleCoordinate(6,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_003_Frac_F_01.mzML", new SampleCoordinate(7,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_005_Frac_F_01.mzML", new SampleCoordinate(8,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_006_Frac_F_01.mzML", new SampleCoordinate(9,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_008_Frac_M_01.mzML", new SampleCoordinate(10,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_009_Frac_M_01.mzML", new SampleCoordinate(11,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_011_Frac_F_01.mzML", new SampleCoordinate(12,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_012_Frac_M_01.mzML", new SampleCoordinate(13,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_013_Frac_F_01.mzML", new SampleCoordinate(14,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_015_Frac_M_01.mzML", new SampleCoordinate(15,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_016_Frac_M_01.mzML", new SampleCoordinate(16,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_018_Frac_M_01.mzML", new SampleCoordinate(17,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_019_Frac_M_01.mzML", new SampleCoordinate(18,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_020_Frac_M_01.mzML", new SampleCoordinate(19,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_021_Frac_F_01.mzML", new SampleCoordinate(20,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_022_Frac_F_01.mzML", new SampleCoordinate(21,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_029_Frac_F_01.mzML", new SampleCoordinate(22,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_030_Frac_M_01.mzML", new SampleCoordinate(23,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_031_Frac_M_01.mzML", new SampleCoordinate(24,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_033_Frac_F_01.mzML", new SampleCoordinate(25,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_038_Frac_M_01.mzML", new SampleCoordinate(26,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_039_Frac_F_01.mzML", new SampleCoordinate(27,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_040_Frac_M_01.mzML", new SampleCoordinate(28,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_045_Frac_M_01.mzML", new SampleCoordinate(29,4));
		sampleKey.put("032922_pancreatitis_1ug_DIA_UPF_047_Frac_F_01.mzML", new SampleCoordinate(30,4));
	}
}
