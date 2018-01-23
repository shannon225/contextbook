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
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TObjectIntHashMap;

public class PerspectiveElibParser {
	public static void main(String[] args) {
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
			library.openFile(new File("/Users/searleb/Documents/school/perspective/normal_combined.elib"));
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
			library.openFile(new File("/Users/searleb/Documents/school/perspective/variable_combined.elib"));
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
			library.openFile(new File("/Users/searleb/Documents/school/perspective/overlap_combined.elib"));
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
