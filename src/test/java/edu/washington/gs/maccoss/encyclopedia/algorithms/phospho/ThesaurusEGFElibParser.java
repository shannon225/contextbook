package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.stat.inference.TestUtils;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideTrie;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SimplePeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.StringUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.BenjaminiHochberg;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedianDouble;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class ThesaurusEGFElibParser {
	public static TreeMap<String, Coordinate> sampleKey=new TreeMap<>();
	
	private static final int numberOfSampleTypes=6;
	private static final int numberOfReplicates=4;
	private static String getSampleName(int i) {
		switch (i) {
		case 1: return "H0+DMSO";
		case 2: return "H0+EGF";
		case 3: return "H4+DMSO";
		case 4: return "H4+EGF";
		case 5: return "H16+DMSO";
		case 6: return "H16+EGF";
		default: return i > 0 && i < 27 ? String.valueOf((char)(i + 64)) : null;
		}
	}
	
	private static final int[] getTrend(float[][] data) {
		int[] r=new int[5];
		for (int i=0; i<r.length; i++) {
			testChange(data, r, i+1);	
		}
		return r;
	}

	private static void testChange(float[][] data, int[] r, int test) {
		double p=TestUtils.pairedTTest(General.toDoubleArray(data[0]), General.toDoubleArray(data[test]));
		if (p<0.001) r[test-1]=((General.mean(data[0])>General.mean(data[test]))?-1:1);
	}
	
	private static final byte H0_VS_H4_DMSO=1;
	private static final byte H4_VS_H16_DMSO=2;
	private static final byte H4_VS_H16_EGF=3;
	private static final byte H4_EGF_VS_DMSO=4;
	private static final byte H16_EGF_VS_DMSO=5;
	private static final byte H0_EGF_VS_DMSO=6;
	private static final byte EGF_VS_DMSO=7;
	
	private static final double getPairedTTest(float[][] data, byte testType) {
		Pair<TDoubleArrayList, TDoubleArrayList> xy=getPairedData(data, testType);	
		return TestUtils.pairedTTest(xy.x.toArray(), xy.y.toArray());
	}
	
	private static final double getFoldChange(float[][] data, byte testType) {
		Pair<TDoubleArrayList, TDoubleArrayList> xy=getPairedData(data, testType);
		
		double[] x=xy.x.toArray();
		double[] y=xy.y.toArray();
		
		if (true) {
			double medianX=QuickMedianDouble.median(x);
			double medianY=QuickMedianDouble.median(y);
			if (medianX==0.0&&medianY==0.0) {
				return Double.NaN;
			} else {
				if (medianX==0.0) {
					return 5;
				} else if (medianY==0.0) {
					return -5;
				} else {
					return Log.log2(General.mean(y)/General.mean(x));
				}
			}
		}
		
		TDoubleArrayList fc=new TDoubleArrayList();
		for (int i=0; i<x.length; i++) {
			if (x[i]!=0.0&&y[i]!=0.0) {
				if (x[i]==0.0) {
					fc.add(100);
				} else if (y[i]==0.0) {
					fc.add(-100);
				} else {
					fc.add(Log.log2(y[i]/x[i]));
				}
			}
		}
		return QuickMedianDouble.median(fc.toArray());
	}
	
	private static final String[] interestingProteins=new String[] {"O14733", "P00533", "P00533", "P00533", "P04049", "P04049", "P04049", "P04626", "P04626", "P04637", "P04637", "P04637", "P05412",
			"P05412", "P05412", "P05412", "P06239", "P07948", "P10398", "P12931", "P12931", "P15056", "P16885", "P16885", "P19174", "P19174", "P21860", "P21860", "P22681", "P23443", "P23443",
			"P27361", "P27361", "P31749", "P31749", "P31749", "P31751", "P31751", "P31751", "P35222", "P35222", "P35222", "P40763", "P40763", "P42224", "P42224", "P42345", "P42345", "P43403",
			"P45983", "P45984", "P45985", "P45985", "P45985", "P49023", "P49840", "P49841", "P49841", "P51692", "P56945", "P56945", "Q00987", "Q01970", "Q01970", "Q02078", "Q02750", "Q02750",
			"Q02750", "Q02750", "Q03135", "Q05397", "Q05397", "Q05397", "Q12778", "Q12778", "Q13153", "Q13164", "Q13480", "Q13480", "Q13480", "Q13555", "Q13557", "Q9NYJ8", "Q9UQC2"};

	private static Pair<TDoubleArrayList, TDoubleArrayList> getPairedData(float[][] data, byte testType) {
		TDoubleArrayList x=new TDoubleArrayList();
		TDoubleArrayList y=new TDoubleArrayList();
		Pair<TDoubleArrayList, TDoubleArrayList> xy=new Pair<TDoubleArrayList, TDoubleArrayList>(x, y);
		
		switch (testType) {
		case H0_VS_H4_DMSO:
			x.addAll(General.toDoubleArray(data[0]));
			y.addAll(General.toDoubleArray(data[2]));
			break;
		case H4_VS_H16_DMSO:
			x.addAll(General.toDoubleArray(data[2]));
			y.addAll(General.toDoubleArray(data[4]));
			break;
		case H0_EGF_VS_DMSO:
			x.addAll(General.toDoubleArray(data[0]));
			y.addAll(General.toDoubleArray(data[1]));
			break;
		case H4_EGF_VS_DMSO:
			x.addAll(General.toDoubleArray(data[2]));
			y.addAll(General.toDoubleArray(data[3]));
			break;
		case H16_EGF_VS_DMSO:
			x.addAll(General.toDoubleArray(data[4]));
			y.addAll(General.toDoubleArray(data[5]));
			break;
		case H4_VS_H16_EGF:
			x.addAll(General.toDoubleArray(data[3]));
			y.addAll(General.toDoubleArray(data[5]));
			break;
		case EGF_VS_DMSO:
			//x.addAll(General.toDoubleArray(data[0]));
			//y.addAll(General.toDoubleArray(data[1]));
			x.addAll(General.toDoubleArray(data[2]));
			y.addAll(General.toDoubleArray(data[3]));
			x.addAll(General.toDoubleArray(data[4]));
			y.addAll(General.toDoubleArray(data[5]));
			break;
		default:
			break;
		}
		return xy;
	}
	
	public static void loadMap() {
		sampleKey.put("2017nov21_bcs_qe1_phospho_1a.mzML", new Coordinate(1, 1));
		sampleKey.put("2017nov21_bcs_qe1_phospho_1b.mzML", new Coordinate(1, 2));
		sampleKey.put("2017nov21_bcs_qe1_phospho_1c.mzML", new Coordinate(1, 3));
		sampleKey.put("2017nov21_bcs_qe1_phospho_1d.mzML", new Coordinate(1, 4));
		sampleKey.put("2017nov21_bcs_qe1_phospho_1e.mzML", new Coordinate(1, 5));
		sampleKey.put("2017nov21_bcs_qe1_phospho_1f.mzML", new Coordinate(1, 6));
		sampleKey.put("2017nov21_bcs_qe1_phospho_2a.mzML", new Coordinate(2, 1));
		sampleKey.put("2017nov21_bcs_qe1_phospho_2b.mzML", new Coordinate(2, 2));
		sampleKey.put("2017nov21_bcs_qe1_phospho_2c.mzML", new Coordinate(2, 3));
		sampleKey.put("2017nov21_bcs_qe1_phospho_2d.mzML", new Coordinate(2, 4));
		sampleKey.put("2017nov21_bcs_qe1_phospho_2e.mzML", new Coordinate(2, 5));
		sampleKey.put("2017nov21_bcs_qe1_phospho_2f.mzML", new Coordinate(2, 6));
		sampleKey.put("2017nov21_bcs_qe1_phospho_3a.mzML", new Coordinate(3, 1));
		sampleKey.put("2017nov21_bcs_qe1_phospho_3b.mzML", new Coordinate(3, 2));
		sampleKey.put("2017nov21_bcs_qe1_phospho_3c.mzML", new Coordinate(3, 3));
		sampleKey.put("2017nov21_bcs_qe1_phospho_3d.mzML", new Coordinate(3, 4));
		sampleKey.put("2017nov21_bcs_qe1_phospho_3e.mzML", new Coordinate(3, 5));
		sampleKey.put("2017nov21_bcs_qe1_phospho_3f.mzML", new Coordinate(3, 6));
		sampleKey.put("2017nov21_bcs_qe1_phospho_4a.mzML", new Coordinate(4, 1));
		sampleKey.put("2017nov21_bcs_qe1_phospho_4b.mzML", new Coordinate(4, 2));
		sampleKey.put("2017nov21_bcs_qe1_phospho_4c.mzML", new Coordinate(4, 3));
		sampleKey.put("2017nov21_bcs_qe1_phospho_4d.mzML", new Coordinate(4, 4));
		sampleKey.put("2017nov21_bcs_qe1_phospho_4e.mzML", new Coordinate(4, 5));
		sampleKey.put("2017nov21_bcs_qe1_phospho_4f.mzML", new Coordinate(4, 6));
	}
	public static final boolean TOTAL_ANALYSIS=true;
	
	public static final boolean MOTIF_ANALYSIS=false;
	public static final boolean ANOVA_ANALYSIS=false;
	public static final boolean HEATMAP_ANALYSIS=true;
	public static final boolean MULTIPLE_FORM_ANALYSIS=false;
	public static final boolean SITE_SPECIFIC_VS_TOTAL_ANALYSIS=false;
	
	public static void main(String[] args) throws Exception {
		
		LibraryFile.OPEN_IN_PLACE=true;
		Logger.PRINT_TO_SCREEN=false;
		loadMap();
		byte targetFoldChangeData=H4_EGF_VS_DMSO;

		PeptideModification mod=PeptideModification.phosphorylation;
		String[] targets=null;//interestingProteins; //null
		
		File[] f=new File("/Volumes/DataBackup/2017nov21_hela_egf").listFiles();
		
		Pair<TreeMap<String,QuantitationLog>, TreeMap<String,QuantitationLog>> quantLogPair=getQuantData(targets, f);
		TreeMap<String, QuantitationLog> totalQuantLog=quantLogPair.x;
		TreeMap<String, QuantitationLog> siteSpecificQuantLog=quantLogPair.y;
		float[][] totalQuantLogNormFactor=getSampleNormalization(totalQuantLog);
		float[][] siteSpecificQuantLogNormFactor=getSampleNormalization(siteSpecificQuantLog);
		
		PeptideMotifTrie motifTrie=new PeptideMotifTrie(totalQuantLog.values(), mod);

		System.out.println("Reading FASTA...");
		//ArrayList<FastaEntryInterface> fasta=FastaReader.readFasta(new File("/Users/searleb/Documents/school/projects/pecandata/UP000005640_9606.fasta"));
		ArrayList<FastaEntryInterface> fasta=FastaReader.readFasta(new File("/Users/searleb/Documents/data/dbs/UP000005640_9606.fasta"));
		motifTrie.addFasta(fasta);
		
		ArrayList<String> siteSpecificPeptides=new ArrayList<>();
		TDoubleArrayList siteSpecificPValues=new TDoubleArrayList();
		TDoubleArrayList siteSpecificAnova=new TDoubleArrayList();
		TDoubleArrayList siteSpecificFC=new TDoubleArrayList();
		for (String peptide : siteSpecificQuantLog.keySet()) {
			QuantitationLog log=siteSpecificQuantLog.get(peptide);
			if (log.getNumMeasurements()<12) continue;
			if (!log.isAtLeastOneCaseFull()) continue;
			
			float[][] data=log.getNormalizedData(siteSpecificQuantLogNormFactor);
			double pValue;
			
			pValue=getPairedTTest(data, EGF_VS_DMSO);
			
			if (Double.isNaN(pValue)) continue;
			if (Double.isInfinite(pValue)) continue;
			if (pValue<0) continue;
			siteSpecificPeptides.add(peptide);
			siteSpecificPValues.add(pValue);
			siteSpecificFC.add(getFoldChange(data, targetFoldChangeData));

			pValue=getANOVAPValue(data);
			siteSpecificAnova.add(pValue);
		}
		double[] siteSpecificAdjustedPValues=BenjaminiHochberg.calculateAdjustedPValues(siteSpecificPValues.toArray());
		double[] siteSpecificAdjustedAnova=BenjaminiHochberg.calculateAdjustedPValues(siteSpecificAnova.toArray());

		ArrayList<String> totalPeptides=new ArrayList<>();
		TDoubleArrayList pValues=new TDoubleArrayList();
		TDoubleArrayList totalAnova=new TDoubleArrayList();
		TDoubleArrayList totalFC=new TDoubleArrayList();
		for (String peptide : totalQuantLog.keySet()) {
			QuantitationLog log=totalQuantLog.get(peptide);
			if (log.getNumMeasurements()<12) continue;
			if (!log.isAtLeastOneCaseFull()) continue;
			
			float[][] data=log.getNormalizedData(totalQuantLogNormFactor);
			double pValue;

			pValue=getPairedTTest(data, EGF_VS_DMSO);
			
			if (Double.isNaN(pValue)) continue;
			if (Double.isInfinite(pValue)) continue;
			if (pValue<0) continue;
			totalPeptides.add(peptide);
			pValues.add(pValue);
			totalFC.add(getFoldChange(data, targetFoldChangeData));

			pValue=getANOVAPValue(data);
			totalAnova.add(pValue);
		}
		double[] totalAdjustedPValues=BenjaminiHochberg.calculateAdjustedPValues(pValues.toArray());
		double[] totalAdjustedAnova=BenjaminiHochberg.calculateAdjustedPValues(totalAnova.toArray());
		
		TreeMap<String, QuantitationLog> primaryQuantLog=TOTAL_ANALYSIS?totalQuantLog:siteSpecificQuantLog;
		float[][] primaryQuantLogNormFactor=TOTAL_ANALYSIS?totalQuantLogNormFactor:siteSpecificQuantLogNormFactor;
		double[] primaryAdjustedPValues=TOTAL_ANALYSIS?totalAdjustedPValues:siteSpecificAdjustedPValues;
		double[] primaryAdjustedAnova=TOTAL_ANALYSIS?totalAdjustedAnova:siteSpecificAdjustedAnova;
		ArrayList<String> primaryPeptides=TOTAL_ANALYSIS?totalPeptides:siteSpecificPeptides;
		
		if (MOTIF_ANALYSIS) {
			TreeSet<String> motifs=new TreeSet<>();

			for (int pep=0; pep<primaryAdjustedPValues.length; pep++) {
				if (primaryAdjustedPValues[pep]<0.05) {
					String peptide=primaryPeptides.get(pep);
					QuantitationLog log=primaryQuantLog.get(peptide);
					if (log.motifMap.size()!=0) {
						motifs.addAll(log.getMotifs());
					}
				}
			}
			for (String motif : motifs) {
				System.out.println(motif);
			}
		}

		if (ANOVA_ANALYSIS) {
			for (int pep=0; pep<primaryAdjustedPValues.length; pep++) {
				String peptide=primaryPeptides.get(pep);
				QuantitationLog log=primaryQuantLog.get(peptide);

				if (primaryAdjustedPValues[pep]<0.05) { // true|| NOTE REPORTS ALL, NOT JUST FDR CORRECTED

					float[][] data=log.getNormalizedData(primaryQuantLogNormFactor);
					double pValue=getANOVAPValue(data);

					System.out.println(log.peptideModSeq+"+"+log.charge+" rt:"+General.mean(log.rtInSecondsList.toArray())+" localized:"+log.isSiteSpecific+" p="+pValue+", FDR="+primaryAdjustedPValues[pep]+" ("+log.toSitesString()+")");
					boolean first=true;
					for (int samp=data.length-1; samp>=0; samp--) {
						if (first) {
							first=false;
						} else {
							System.out.print('\t');
						}
						System.out.print(getSampleName(samp+1));
					}
					System.out.println();
					for (int rep=0; rep<data[0].length; rep++) {
						first=true;
						for (int samp=data.length-1; samp>=0; samp--) {
							if (first) {
								first=false;
							} else {
								System.out.print('\t');
							}
							if (data[samp][rep]>0) {
								System.out.print(data[samp][rep]);
							} else {
								System.out.print(0);
							}
						}
						System.out.println();
					}
					System.out.println();
				}
			}
		}

		if (HEATMAP_ANALYSIS) {
			System.out.println();
			System.out.print("Peptide\tProtein\tlocalized\tANOVA\tPAIRED\tAKT\tLAKT\tMTOR\tMAPK");
			for (int i=0; i<6; i++) {
				System.out.print('\t');
				System.out.print(getSampleName(i+1));
			}
			System.out.println();

			ArrayList<String> flagged=new ArrayList<>();
			for (int pep=0; pep<primaryAdjustedPValues.length; pep++) {
				if (primaryAdjustedPValues[pep]<0.05) { // REMOVE TRUE
					String peptide=primaryPeptides.get(pep);
					QuantitationLog log=primaryQuantLog.get(peptide);

					float[][] data=log.getNormalizedData(primaryQuantLogNormFactor);

					//double pValue=getANOVAPValue(data);
					//pValue=General.min(new double[] {getPairedTTest(data, H0_EGF_VS_DMSO), getPairedTTest(data, H4_EGF_VS_DMSO), getPairedTTest(data, H16_EGF_VS_DMSO)});
					System.out.print(log.peptideModSeq+"\t"+log.toSitesString()+"\t"+log.isSiteSpecific+"\t"+primaryAdjustedAnova[pep]+"\t"+totalAdjustedPValues[pep]);
					System.out.print("\t"+(log.doesMotifMatch("R.R..[ST].....")));
					System.out.print("\t"+(log.doesMotifMatch("..R..[ST].....")));
					System.out.print("\t"+(log.doesMotifMatch(".....[ST][FLW]....")));
					System.out.print("\t"+(log.doesMotifMatch(".....[ST]P....")));
					for (int i=0; i<data.length; i++) {
						System.out.print('\t');
						float[] replicates=data[i];
						System.out.print(getMedianValue(replicates));
					}
					System.out.println();

					//if (primaryAdjustedPValues[pep]<pValue) {
					//	flagged.add(peptide);
					//}
				}
			}
			
			if (flagged.size()>0) {
				System.out.println("FLAGGED!");
				for (String string : flagged) {
					System.out.println("\t"+string);
				}
			}
		}
		
		if (MULTIPLE_FORM_ANALYSIS) {
			TreeMap<String, TDoubleArrayList> pvalueMap=new TreeMap<>();
			TreeMap<String, TFloatArrayList> rtMap=new TreeMap<>();
			double[] values=siteSpecificFC.toArray(); //adjustedPValues;
			for (int pep=0; pep<values.length; pep++) {
				if (Double.isNaN(values[pep])) continue;
				
				String peptide=siteSpecificPeptides.get(pep);
				QuantitationLog log=primaryQuantLog.get(peptide);
				String key=PeptideUtils.getPeptideSeq(peptide);
				
				TDoubleArrayList list=pvalueMap.get(key);
				TFloatArrayList rtList=rtMap.get(key);
				if (list==null) {
					list=new TDoubleArrayList();
					pvalueMap.put(key, list);
					rtList=new TFloatArrayList();
					rtMap.put(key, rtList);
				}
				list.add(values[pep]);
				rtList.add(General.mean(log.rtInSecondsList.toArray()));
				
			}
			
			for (Entry<String, TDoubleArrayList> entry : pvalueMap.entrySet()) {
				if (entry.getValue().size()>1) {
					double minPValue=Double.MAX_VALUE;
					double maxPValue=-Double.MAX_VALUE;
					double minRT=0;
					double maxRT=0;
					double[] pvalues=entry.getValue().toArray();
					float[] rts=rtMap.get(entry.getKey()).toArray();
					for (int j=0; j<pvalues.length; j++) {
						if (pvalues[j]>maxPValue) {
							maxPValue=pvalues[j];
							maxRT=rts[j];
						}
						if (pvalues[j]<minPValue) {
							minPValue=pvalues[j];
							minRT=rts[j];
						}
					}
					System.out.println(entry.getKey()+"\t"+entry.getValue().size()+"\t"+minPValue+"\t"+maxPValue+"\t"+(maxRT-minRT));
				}
			}
			System.out.println(pvalueMap.size()+" Total forms");
		}
		
		if (SITE_SPECIFIC_VS_TOTAL_ANALYSIS) {
			TreeMap<String, double[]> pvalueMap=new TreeMap<>();
			for (int pep=0; pep<totalAdjustedPValues.length; pep++) {
				String peptide=primaryPeptides.get(pep);
				double[] list=pvalueMap.get(peptide);
				if (list==null) {
					list=new double[2];
					Arrays.fill(list, -1);
					pvalueMap.put(peptide, list);
				}
				list[0]=totalAdjustedPValues[pep];
			}

			for (int pep=0; pep<siteSpecificAdjustedPValues.length; pep++) {
				String peptide=siteSpecificPeptides.get(pep);
				double[] list=pvalueMap.get(peptide);
				if (list==null) {
					list=new double[2];
					Arrays.fill(list, -1);
					pvalueMap.put(peptide, list);
				}
				list[1]=siteSpecificAdjustedPValues[pep];
			}
			
			int count=0;
			for (Entry<String, double[]> entry : pvalueMap.entrySet()) {
				double[] pair=entry.getValue();
				if (pair[0]>=0.25&&pair[1]>=0&&pair[1]<0.01) {
					count++;
					QuantitationLog log=siteSpecificQuantLog.get(entry.getKey());
					System.out.println(entry.getKey()+"\t"+log.charge+"\t"+General.mean(log.rtInSecondsList.toArray())+"\t"+General.max(log.localizationScores.toArray())+"\t"+General.min(log.localizationScores.toArray())+"\t"+pair[0]+"\t"+pair[1]);
				}
			}
			System.out.println(count+" Total forms");
		}
	}

	public static float getMedianValue(float[] replicates) {
		return QuickMedian.median(replicates.clone());
		
//		float[] clone=replicates.clone();
//		Arrays.sort(clone);
//		int index;
//		if (clone.length%2==0) {
//			// even     0,1,2,3: floor(4/2)=2 (also need 1)
//			// even 0,1,2,3,4,5: floor(6/2)=3 (also need 2)
//			index=(clone.length/2);
//			if (clone[index-1]>0) {
//				return (clone[index-1]+clone[index])/2.0f;
//			} else if (clone[index]>0) {
//				return clone[index];
//			}
//		} else {
//			//     0,1,2: floor(3/2)=1
//			// 0,1,2,3,4: floor(5/2)=2
//			index=(clone.length/2);
//			if (clone[index]>0) return clone[index];
//		}
//		for (index=index+1; index<clone.length; index++) {
//			if (clone[index]>0) return clone[index];
//		}
//		return 0;
	}

	private static Pair<TreeMap<String, QuantitationLog>, TreeMap<String, QuantitationLog>> getQuantData(String[] targets, File[] f) throws IOException, SQLException {
		TreeMap<String, QuantitationLog> quantLog=new TreeMap<>();
		TreeMap<String, QuantitationLog> siteSpecificQuantLog=new TreeMap<>();
		for (File file : f) {
			if (file.getName().endsWith(CASiLJobData.THESAURUS_REPORT_FILE_SUFFIX)&&!file.getName().startsWith(".")) {
				System.out.println("Parsing "+file.getName()+"...");
				LibraryFile library=new LibraryFile();
				library.openFile(file);
				

				Connection c=library.getConnection();
				Statement s=c.createStatement();
				ResultSet rs = s.executeQuery("select pep.PrecursorCharge, pep.PeptideModSeq, pep.PeptideSeq, pep.SourceFile, max(pep.LocalizedIntensity), max(pep.TotalIntensity), pep.IsSiteSpecific, pep.RTInSecondsCenter,pep.localizationScore,"+
						"group_concat(p.ProteinAccession, '" + PSMData.ACCESSION_TOKEN + "') as ProteinAccessions " +
						"from " +
						"peptidelocalizations pep " +
						"left join peptidetoprotein p " +
						"where " +
						"pep.PeptideSeq = p.PeptideSeq " +
						"group by pep.rowid;"
				);
				while (rs.next()) {
					byte precursorCharge=(byte)rs.getInt(1);
					String peptideModSeq=rs.getString(2);
					String peptideSeq=rs.getString(3);
					String sourceFile=rs.getString(4);
					float localizedIntensity=rs.getFloat(5);
					float totalIntensity=rs.getFloat(6);
					boolean isSiteSpecific=rs.getBoolean(7);
					float rtInSeconds=rs.getFloat(8);
					float localizationScore=rs.getFloat(9);
					String proteinToken=rs.getString(10);
					//HashSet<String> accessions=PSMData.stringToAccessions(proteinToken);
					
					boolean keeper=false;
					if (targets==null) {
						keeper=true;
					} else {
						for (int i=0; i<targets.length; i++) {
							if (proteinToken.indexOf(targets[i])>=0) {
								keeper=true;
								break;
							}
						}
					}
					
					if (keeper) {
						Coordinate coord=sampleKey.get(sourceFile);
						if (coord==null) {
							System.out.println("FAILED TO FIND SAMPLE: "+sourceFile);
							System.exit(1);
						}

						QuantitationLog log=quantLog.get(peptideModSeq);
						QuantitationLog siteLog=siteSpecificQuantLog.get(peptideModSeq);
						if (log==null) {
							log=new QuantitationLog(proteinToken, peptideModSeq, precursorCharge);
							quantLog.put(peptideModSeq, log);
							siteLog=new QuantitationLog(proteinToken, peptideModSeq, precursorCharge);
							siteSpecificQuantLog.put(peptideModSeq, siteLog);
						}
						log.addIntensity(coord, totalIntensity, rtInSeconds, localizationScore, isSiteSpecific);
						siteLog.addIntensity(coord, localizedIntensity, rtInSeconds, localizationScore, isSiteSpecific);
						
					}
				}
				rs.close();
				s.close();
				c.close();
			}
		}
		return new Pair<TreeMap<String,QuantitationLog>, TreeMap<String,QuantitationLog>>(quantLog, siteSpecificQuantLog);
	}
	
	private static float[][] getSampleNormalization(TreeMap<String,QuantitationLog> logs) {
		float[][] totals=new float[numberOfSampleTypes][];
		for (int i=0; i<totals.length; i++) {
			totals[i]=new float[numberOfReplicates];
		}
		for (QuantitationLog log : logs.values()) {
			float[][] data=log.getData();
			for (int i=0; i<data.length; i++) {
				totals[i]=General.add(totals[i], data[i]);
			}
		}
		float averageTotal=General.mean(totals);
		for (int i=0; i<totals.length; i++) {
			for (int j=0; j<totals[i].length; j++) {
				totals[i][j]=averageTotal/totals[i][j];
			}
		}
		return totals;
	}
	
	private static double getANOVAPValue(float[][] data) {
		ArrayList<double[]> classes=new ArrayList<>();
		for (int i=0; i<data.length; i++) {
			classes.add(General.toDoubleArray(data[i]));
		}
		return TestUtils.oneWayAnovaPValue(classes);
	}
	
	public static class PeptideMotifTrie extends PeptideTrie<QuantitationLog> {
		PeptideModification mod;
		public PeptideMotifTrie(Collection<QuantitationLog> entries, PeptideModification mod) {
			super(entries);
			this.mod=mod;
		}

		@Override
		protected void processMatch(FastaEntryInterface fasta, QuantitationLog entry, int start) {
			int[] indicies=PeptideUtils.getModIndicies(entry.getPeptideModSeq(), mod.getNominalMass());
			for (int i=0; i<indicies.length; i++) {
				int beginIndex=start+indicies[i]-6;
				int leftPad=Math.max(0, -beginIndex);
				int endIndex=start+indicies[i]+5;
				int rightPad=Math.max(0, endIndex-fasta.getSequence().length());
				String motif=(StringUtils.getPad(leftPad, 'X'))+(fasta.getSequence().substring(beginIndex+leftPad, endIndex-rightPad))+(StringUtils.getPad(rightPad, 'X'));
				entry.addMotif(fasta.getAccession(), start+indicies[i], motif);
			}
		}
	}
	
	public static class QuantitationLog extends SimplePeptidePrecursor {
		boolean isSiteSpecific=false;
		final byte charge;
		TreeMap<String, TIntObjectHashMap<String>> motifMap=new TreeMap<String, TIntObjectHashMap<String>>();
		final String protein;
		final String peptideModSeq;
		final TObjectFloatHashMap<Coordinate> intensities=new TObjectFloatHashMap<>();
		final TFloatArrayList rtInSecondsList=new TFloatArrayList();
		final TFloatArrayList localizationScores=new TFloatArrayList();
		
		public QuantitationLog(String protein, String peptideModSeq, byte charge) {
			super(peptideModSeq, charge);
			this.protein=protein;
			this.peptideModSeq=peptideModSeq;
			this.charge=charge;
		}

		public String toSitesString() {
			StringBuilder sb=new StringBuilder();
			for (Entry<String, TIntObjectHashMap<String>> entry : motifMap.entrySet()) {
				if (sb.length()>0) sb.append(";");
				sb.append(entry.getKey()+"("+General.toString(entry.getValue().keys())+")");
			}
			if (sb.length()==0) return protein+"(?)";
			return sb.toString();
		}
		
		public HashSet<String> getMotifs() {
			HashSet<String> set=new HashSet<String>();
			for (TIntObjectHashMap<String> map : motifMap.values()) {
				set.addAll(map.valueCollection());
			}
			return set;
		}
		
		public void addMotif(String accession, int index, String motif) {
			TIntObjectHashMap<String> map=motifMap.get(accession);
			if (map==null) {
				map=new TIntObjectHashMap<>();
				motifMap.put(accession, map);
			}
			map.put(index, motif);
		}
		
		public boolean doesMotifMatch(String s) {
			for (TIntObjectHashMap<String> map : motifMap.values()) {
				for (String string : map.valueCollection()) {
					if (string.matches(s))
						return true;
				}
			}
			return false;
		}

		public int getNumMeasurements() {
			return intensities.size();
		}
		
		public boolean isAtLeastOneCaseFull() {
			float[][] data=getData();
			for (int i=0; i<data.length; i++) {
				boolean full=true;
				for (int j=0; j<data[i].length; j++) {
					if (data[i][j]==0.0f) {
						full=false;
						break;
					}
				}
				if (full) return true;
			}
			return false;
		}
		
		public void addIntensity(Coordinate c, float intensity, float rtInSeconds, float localizationScore, boolean isSiteSpecific) {
			intensities.adjustOrPutValue(c, intensity, intensity);
			if (isSiteSpecific) {
				this.isSiteSpecific=true;
			}
			this.rtInSecondsList.add(rtInSeconds);
			this.localizationScores.add(localizationScore);
		}
		public float[][] getData() {
			float[][] results=new float[numberOfSampleTypes][];
			for (int i=0; i<results.length; i++) {
				results[i]=new float[numberOfReplicates];
			}
			intensities.forEachEntry(new TObjectFloatProcedure<Coordinate>() {
				@Override
				public boolean execute(Coordinate a, float b) {
					results[a.sample-1][a.replicate-1]+=b;
					return true;
				}
			});
			
			return results;
		}
		public float[][] getNormalizedData(float[][] normalizationMatrix) {
			float[][] data=getData();
			for (int i=0; i<data.length; i++) {
				data[i]=General.multiply(data[i], normalizationMatrix[i]);
			}
			
			float[][] normalized=new float[data.length][];
			for (int i=0; i<normalized.length; i++) {
				normalized[i]=new float[data[i].length];
			}

			float grandTotal=0.0f;
			for (int rep=0; rep<data[0].length; rep++) {
				float total=0.0f;
				for (int samp=0; samp<data.length; samp++) {
					total+=data[samp][rep];
				}
				for (int samp=0; samp<data.length; samp++) {
					normalized[samp][rep]=data[samp][rep]/total;
				}
				grandTotal+=total;
			}
			grandTotal=grandTotal/data[0].length;
			for (int rep=0; rep<data[0].length; rep++) {
				for (int samp=0; samp<data.length; samp++) {
					normalized[samp][rep]=normalized[samp][rep]*grandTotal;
				}
			}
			
			return normalized;
		}
	}
	
	public static class Coordinate {
		final int replicate;
		final int sample;
		
		public Coordinate(int replicate, int sample) {
			this.replicate=replicate;
			this.sample=sample;
		}
		
		@Override
		public boolean equals(Object obj) {
			if (obj==null) return false;
			if (!(obj instanceof Coordinate)) return false;
			return hashCode()==obj.hashCode();
		}
		
		@Override
		public int hashCode() {
			return sample+replicate*1000;
		}
	}
}
