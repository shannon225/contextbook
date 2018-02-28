package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TObjectIntHashMap;

public class PerspectiveElibParser {
	public static void main(String[] args) {
		HashMap<String, Pair<TFloatArrayList, TFloatArrayList>> data=getPRM();
		HashMap<String, Pair<TFloatArrayList, TFloatArrayList>> ddaData=getDDA();
		HashMap<String, TFloatArrayList[]> normalDiaData=getDIA(reportNormal);
		HashMap<String, TFloatArrayList[]> variableDiaData=getDIA(reportVariable);
		HashMap<String, TFloatArrayList[]> overlapDiaData=getDIA(reportOverlap);

		String[] keys=data.keySet().toArray(new String[data.size()]);
		Arrays.sort(keys);
		System.out.println("peptide,intensity,prm_ratio,cv,log intensity,delta_dda,delta_normal_dia,delta_variable_dia,delta_overlap_dia");
		for (String key : keys) {
			Pair<TFloatArrayList, TFloatArrayList> prmMeasurements=data.get(key);
			float[] xs=prmMeasurements.x.toArray().clone();
			float x=QuickMedian.median(xs);
			float[] ys=prmMeasurements.y.toArray().clone();
			float y=QuickMedian.median(ys);
			float prmRatio=protectedLog2(y/x);

			float[] lmA=General.divide(xs, x);
			float[] lmB=General.divide(ys, y);
			float cv=General.stdev(General.concatenate(lmA, lmB));

			Pair<TFloatArrayList, TFloatArrayList> ddaMeasurements=ddaData.get(key);
			Float ddaRatio=ddaMeasurements==null?null:protectedLog2(QuickMedian.median(ddaMeasurements.y.toArray().clone())/QuickMedian.median(ddaMeasurements.x.toArray().clone()));

			TFloatArrayList[] diaMeasurements=normalDiaData.get(key);
			Float normalDIARatio=diaMeasurements==null?null:protectedLog2(QuickMedian.median(diaMeasurements[1].toArray().clone())/QuickMedian.median(diaMeasurements[0].toArray().clone()));

			diaMeasurements=variableDiaData.get(key);
			Float variableDIARatio=diaMeasurements==null?null:protectedLog2(QuickMedian.median(diaMeasurements[1].toArray().clone())/QuickMedian.median(diaMeasurements[0].toArray().clone()));
			
			diaMeasurements=overlapDiaData.get(key);
			Float overlapDIARatio=diaMeasurements==null?null:protectedLog2(QuickMedian.median(diaMeasurements[1].toArray().clone())/QuickMedian.median(diaMeasurements[0].toArray().clone()));
			
			System.out.print(key);
			System.out.print(',');
			System.out.print(Math.max(x,y));
			System.out.print(',');
			System.out.print(prmRatio);
			System.out.print(',');
			System.out.print(cv);
			System.out.print(',');
			System.out.print(Log.log10(Math.max(x,y)));
			System.out.print(',');
			System.out.print(ddaRatio==null?"":(ddaRatio-prmRatio));
			System.out.print(',');
			System.out.print(normalDIARatio==null?"":(normalDIARatio-prmRatio));
			System.out.print(',');
			System.out.print(variableDIARatio==null?"":(variableDIARatio-prmRatio));
			System.out.print(',');
			System.out.print(overlapDIARatio==null?"":(overlapDIARatio-prmRatio));
			System.out.println();
		}
	}
	
	private static Float protectedLog2(float v) {
		if (Float.isNaN(v)) return 0.0f;
		if (Float.isInfinite(v)) {
			if (v>0) {
				return 8.0f;
			} else {
				return -8.0f;
			}
		}
		if (v<=0) return -8.0f;
		
		return Log.log2(v);
	}
	
	public static void oldmain(String[] args) {
		HashMap<String, Pair<TFloatArrayList, TFloatArrayList>> data=getPRM();
		HashMap<String, TFloatArrayList[]> diaData=getDIA();
		
		String[] keys=data.keySet().toArray(new String[data.size()]);
		Arrays.sort(keys);
		for (String key : keys) {
			Pair<TFloatArrayList, TFloatArrayList> measurements=data.get(key);
			float[] a=measurements.x.toArray();
			float[] b=measurements.y.toArray();
			
			float meanA=General.mean(a);
			float[] lmA=General.divide(a, meanA);
			float meanB=General.mean(b);
			float[] lmB=General.divide(b, meanB);
			float cv=General.stdev(General.concatenate(lmA, lmB));
			float targetRatio=Log.log2(meanB/meanA);
			
			switch (mode) {
				case reportPRM:
					//System.out.println(key+"\t"+targetRatio);//+"\t"+cv+"\t"+Math.max(meanA,  meanB));
					System.out.println(key+"\t"+Math.max(meanA,  meanB));
					break;

				default:
					// DIA
					TFloatArrayList[] values=diaData.get(key);
					if (values==null) {
						System.out.println(key);
						break;
					}
					float[] diaA=values[0].toArray();
					float[] diaB=values[1].toArray();
					float diameanA=Math.max(1, General.mean(diaA));
					float[] dialmA=General.divide(diaA, diameanA);
					float diameanB=Math.max(1, General.mean(diaB));
					float[] dialmB=General.divide(diaB, diameanB);
					float diacv=General.stdev(General.concatenate(dialmA, dialmB));
					float diatargetRatio=Log.log2(diameanB/diameanA);

					//System.out.println(key+"\t"+(targetRatio-diatargetRatio)+"\t"+diacv+"\t"+Math.max(diameanA,  diameanB));
					System.out.println(key+"\t"+targetRatio+"\t"+diatargetRatio);
					break;
			}
		}
	}
	
	private static final int mode=0;
	private static final int reportPRM=0;
	private static final int reportNormal=1;
	private static final int reportVariable=2;
	private static final int reportOverlap=3;

	private static HashMap<String, TFloatArrayList[]> getDIA() {
		return getDIA(mode);
	}

	private static HashMap<String, TFloatArrayList[]> getDIA(int mode) {
		HashMap<String, TFloatArrayList[]> diaData=null;
		switch (mode) {
			case reportNormal:
				diaData=loadNormalDIA();
				break;
			case reportVariable:
				diaData=loadVariableDIA();
				break;
			case reportOverlap:
				diaData=loadOverlapDIA();
				break;

			default:
				break;
		}
		return diaData;
	}

	private static HashMap<String, Pair<TFloatArrayList, TFloatArrayList>> getDDA() {
		AminoAcidConstants constants=new AminoAcidConstants();
		HashMap<String, Pair<TFloatArrayList, TFloatArrayList>> data=new HashMap<>(); 
		TableParser.parseCSV(new File("/Users/searleb/Documents/school/perspective/DDA peak areas.csv"), new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String peptideModSeq=row.get("peptidemodseq");
				peptideModSeq=PeptideUtils.getCorrectedMasses(peptideModSeq, constants);
				
				String sample=row.get("sample");
				float intensity=1.0f+Float.parseFloat(row.get("intensity")); // one pseudo count
				Pair<TFloatArrayList, TFloatArrayList> measurements=data.get(peptideModSeq);
				if (measurements==null) {
					measurements=new Pair<TFloatArrayList, TFloatArrayList>(new TFloatArrayList(), new TFloatArrayList());
					data.put(peptideModSeq, measurements);
				}
				
				if ("6b".equals(sample)) {
					measurements.x.add(intensity);
				} else {
					measurements.y.add(intensity);
				}
			}
		});
		return data;
	}

	private static HashMap<String, Pair<TFloatArrayList, TFloatArrayList>> getPRM() {
		AminoAcidConstants constants=new AminoAcidConstants();
		HashMap<String, Pair<TFloatArrayList, TFloatArrayList>> data=new HashMap<>(); 
		TableParser.parseCSV(new File("/Users/searleb/Documents/school/perspective/PRM Peak Areas.csv"), new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String peptideModSeq=row.get("Peptide Modified Sequence Monoisotopic Masses");
				peptideModSeq=PeptideUtils.getCorrectedMasses(peptideModSeq, constants);
				
				String sample=row.get("sample");
				float intensity=1.0f+Float.parseFloat(row.get("Total Area Fragment")); // one pseudo count
				Pair<TFloatArrayList, TFloatArrayList> measurements=data.get(peptideModSeq);
				if (measurements==null) {
					measurements=new Pair<TFloatArrayList, TFloatArrayList>(new TFloatArrayList(), new TFloatArrayList());
					data.put(peptideModSeq, measurements);
				}
				
				if ("6b".equals(sample)) {
					measurements.x.add(intensity);
				} else {
					measurements.y.add(intensity);
				}
			}
		});
		return data;
	}
	
	public static HashMap<String, TFloatArrayList[]> loadNormalDIA() {
		try {
			LibraryFile library=new LibraryFile();
			library.openFile(new File("/Users/searleb/Documents/school/perspective/dda_library/normal_combined.elib"));
			TObjectIntHashMap<String> indexBySampleNames=new TObjectIntHashMap<>();
			indexBySampleNames.put("2017dec27_normal_dia_6b_rep1.mzML", 0);
			indexBySampleNames.put("2017dec27_normal_dia_6b_rep2.mzML", 0);
			indexBySampleNames.put("2017dec27_normal_dia_6b_rep3.mzML", 0);
			indexBySampleNames.put("2017dec27_normal_dia_6e_rep1.mzML", 1);
			indexBySampleNames.put("2017dec27_normal_dia_6e_rep2.mzML", 1);
			indexBySampleNames.put("2017dec27_normal_dia_6e_rep3.mzML", 1);
			return loadMap(library, indexBySampleNames);
		} catch (Exception e) {
			e.printStackTrace();
			throw new EncyclopediaException(e.getMessage());
		}
	}
	
	public static HashMap<String, TFloatArrayList[]> loadVariableDIA() {
		try {
			LibraryFile library=new LibraryFile();
			library.openFile(new File("/Users/searleb/Documents/school/perspective/dda_library/variable_combined.elib"));
			TObjectIntHashMap<String> indexBySampleNames=new TObjectIntHashMap<>();
			indexBySampleNames.put("2017dec27_variable_dia_6b_rep1.mzML", 0);
			indexBySampleNames.put("2017dec27_variable_dia_6b_rep2.mzML", 0);
			indexBySampleNames.put("2017dec27_variable_dia_6b_rep3.mzML", 0);
			indexBySampleNames.put("2017dec27_variable_dia_6e_rep1.mzML", 1);
			indexBySampleNames.put("2017dec27_variable_dia_6e_rep2.mzML", 1);
			indexBySampleNames.put("2017dec27_variable_dia_6e_rep3.mzML", 1);
			return loadMap(library, indexBySampleNames);
		} catch (Exception e) {
			e.printStackTrace();
			throw new EncyclopediaException(e.getMessage());
		}
	}
	
	public static HashMap<String, TFloatArrayList[]> loadOverlapDIA() {
		try {
			LibraryFile library=new LibraryFile();
			library.openFile(new File("/Users/searleb/Documents/school/perspective/dda_library/overlap_combined.elib"));
			TObjectIntHashMap<String> indexBySampleNames=new TObjectIntHashMap<>();
			indexBySampleNames.put("2017dec27_overlap_dia_6b_rep1.mzML", 0);
			indexBySampleNames.put("2017dec27_overlap_dia_6b_rep2.mzML", 0);
			indexBySampleNames.put("2017dec27_overlap_dia_6b_rep3.mzML", 0);
			indexBySampleNames.put("2017dec27_overlap_dia_6e_rep1.mzML", 1);
			indexBySampleNames.put("2017dec27_overlap_dia_6e_rep2.mzML", 1);
			indexBySampleNames.put("2017dec27_overlap_dia_6e_rep3.mzML", 1);
			return loadMap(library, indexBySampleNames);
		} catch (Exception e) {
			e.printStackTrace();
			throw new EncyclopediaException(e.getMessage());
		}
	}
	
	public static HashMap<String, TFloatArrayList[]> loadMap(LibraryFile library, TObjectIntHashMap<String> indexBySampleNames) throws Exception {
		int numGroups=General.max(indexBySampleNames.values())+1;
		HashMap<String, TFloatArrayList[]> data=new HashMap<>(); 
		
		Connection c=library.getConnection();
		Statement s=c.createStatement();
		
		ResultSet rs = s.executeQuery("select " +
				"pep.PeptideModSeq, " +
				"pep.SourceFile, " +
				"pep.TotalIntensity " +
				"from " +
				"peptidequants pep "
		);
		
		while (rs.next()) {
			String peptideModSeq=rs.getString(1);
			String sourceFile=rs.getString(2);
			float totalIntensity=rs.getFloat(3);
			
			TFloatArrayList[] vals=data.get(peptideModSeq);
			if (vals==null) {
				vals=new TFloatArrayList[numGroups];
				for (int i=0; i<vals.length; i++) {
					vals[i]=new TFloatArrayList();
				}
				data.put(peptideModSeq, vals);
			}
			
			int index=indexBySampleNames.get(sourceFile);
			vals[index].add(totalIntensity);
		}
		
		s.close();
		c.close();
		
		return data;
	}
}
