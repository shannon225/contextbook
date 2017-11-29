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
import java.util.HashMap;
import java.util.Map.Entry;
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
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class ThesaurusEGFElibParser {
	public static HashMap<String, Coordinate> sampleKey=new HashMap<>();
	
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
		default:
			break;
		}
		return xy;
	}
	
	public static void loadMap() {
		// FIXME
		sampleKey.put("22jun2016_mcf7_phospho_1a.mzML", new Coordinate(1, 1));
		sampleKey.put("22jun2016_mcf7_phospho_1b.mzML", new Coordinate(1, 2));
		sampleKey.put("22jun2016_mcf7_phospho_1c.mzML", new Coordinate(1, 3));
		sampleKey.put("22jun2016_mcf7_phospho_1d.mzML", new Coordinate(1, 4));
		sampleKey.put("22jun2016_mcf7_phospho_1e.mzML", new Coordinate(1, 5));
		sampleKey.put("22jun2016_mcf7_phospho_1f.mzML", new Coordinate(1, 6));
		sampleKey.put("22jun2016_mcf7_phospho_2a.mzML", new Coordinate(2, 1));
		sampleKey.put("22jun2016_mcf7_phospho_2b.mzML", new Coordinate(2, 2));
		sampleKey.put("22jun2016_mcf7_phospho_2c.mzML", new Coordinate(2, 3));
		sampleKey.put("22jun2016_mcf7_phospho_2d.mzML", new Coordinate(2, 4));
		sampleKey.put("22jun2016_mcf7_phospho_2e.mzML", new Coordinate(2, 5));
		sampleKey.put("22jun2016_mcf7_phospho_2f.mzML", new Coordinate(2, 6));
		sampleKey.put("22jun2016_mcf7_phospho_3a_160627233451.mzML", new Coordinate(3, 1));
		sampleKey.put("22jun2016_mcf7_phospho_3b_160627142134.mzML", new Coordinate(3, 2));
		sampleKey.put("22jun2016_mcf7_phospho_3c_160628015316.mzML", new Coordinate(3, 3));
		sampleKey.put("22jun2016_mcf7_phospho_3d_160627211625.mzML", new Coordinate(3, 4));
		sampleKey.put("22jun2016_mcf7_phospho_3e_160627185757.mzML", new Coordinate(3, 5));
		sampleKey.put("22jun2016_mcf7_phospho_3f_160627163930.mzML", new Coordinate(3, 6));
		sampleKey.put("22jun2016_mcf7_phospho_4a_160627082406.mzML", new Coordinate(4, 1));
		sampleKey.put("22jun2016_mcf7_phospho_4b_160627034715.mzML", new Coordinate(4, 2));
		sampleKey.put("22jun2016_mcf7_phospho_4c_160627060541.mzML", new Coordinate(4, 3));
		sampleKey.put("22jun2016_mcf7_phospho_4d_160626205159.mzML", new Coordinate(4, 4));
		sampleKey.put("22jun2016_mcf7_phospho_4e_160626231025.mzML", new Coordinate(4, 5));
		sampleKey.put("22jun2016_mcf7_phospho_4f_160627012850.mzML", new Coordinate(4, 6));
	}
	public static final boolean TOTAL_ANALYSIS=false;
	
	public static final boolean MOTIF_ANALYSIS=false;
	public static final boolean ANOVA_ANALYSIS=false;
	public static final boolean HEATMAP_ANALYSIS=false;
	public static final boolean MULTIPLE_FORM_ANALYSIS=true;
	public static final boolean SITE_SPECIFIC_VS_TOTAL_ANALYSIS=false;
	
	public static void main(String[] args) throws Exception {
		
		LibraryFile.OPEN_IN_PLACE=true;
		Logger.PRINT_TO_SCREEN=false;
		loadMap();
		byte targetFoldChangeData=H4_EGF_VS_DMSO;

		PeptideModification mod=PeptideModification.phosphorylation;
		String[] targets=null;
		
		File[] f=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs").listFiles(); // FIXME
		
		Pair<HashMap<String,QuantitationLog>, HashMap<String,QuantitationLog>> quantLogPair=getQuantData(targets, f);
		HashMap<String, QuantitationLog> totalQuantLog=quantLogPair.x;
		HashMap<String, QuantitationLog> siteSpecificQuantLog=quantLogPair.y;
		
		HashMap<String, QuantitationLog> primaryQuantLog=TOTAL_ANALYSIS?totalQuantLog:siteSpecificQuantLog;
		
		PeptideMotifTrie motifTrie=new PeptideMotifTrie(primaryQuantLog.values(), mod);

		System.out.println("Reading FASTA...");
		ArrayList<FastaEntryInterface> fasta=FastaReader.readFasta(new File("/Users/searleb/Documents/school/projects/pecandata/UP000005640_9606.fasta"));
		motifTrie.addFasta(fasta);
		
		ArrayList<String> siteSpecificPeptides=new ArrayList<>();
		TDoubleArrayList siteSpecificPValues=new TDoubleArrayList();
		TDoubleArrayList siteSpecificFC=new TDoubleArrayList();
		for (String peptide : siteSpecificQuantLog.keySet()) {
			QuantitationLog log=siteSpecificQuantLog.get(peptide);
			if (log.getNumMeasurements()<12) continue;
			if (!log.isAtLeastOneCaseFull()) continue;
			
			float[][] data=log.getNormalizedData();
			double pValue;
			
			pValue=getANOVAPValue(data);
			
			if (Double.isNaN(pValue)) continue;
			if (Double.isInfinite(pValue)) continue;
			if (pValue<0) continue;
			siteSpecificPeptides.add(peptide);
			siteSpecificPValues.add(pValue);
			siteSpecificFC.add(getFoldChange(data, targetFoldChangeData));
		}
		double[] siteSpecificAdjustedPValues=BenjaminiHochberg.calculateAdjustedPValues(siteSpecificPValues.toArray());

		ArrayList<String> peptides=new ArrayList<>();
		TDoubleArrayList pValues=new TDoubleArrayList();
		TDoubleArrayList totalFC=new TDoubleArrayList();
		for (String peptide : totalQuantLog.keySet()) {
			QuantitationLog log=totalQuantLog.get(peptide);
			if (log.getNumMeasurements()<12) continue;
			if (!log.isAtLeastOneCaseFull()) continue;
			
			float[][] data=log.getNormalizedData();
			double pValue;
			
			pValue=getANOVAPValue(data);
			
			if (Double.isNaN(pValue)) continue;
			if (Double.isInfinite(pValue)) continue;
			if (pValue<0) continue;
			peptides.add(peptide);
			pValues.add(pValue);
			totalFC.add(getFoldChange(data, targetFoldChangeData));
		}
		double[] totalAdjustedPValues=BenjaminiHochberg.calculateAdjustedPValues(pValues.toArray());
		
		if (MOTIF_ANALYSIS) {
			TreeSet<String> motifs=new TreeSet<>();

			for (int pep=0; pep<totalAdjustedPValues.length; pep++) {
				if (totalAdjustedPValues[pep]<0.05) {
					String peptide=peptides.get(pep);
					QuantitationLog log=primaryQuantLog.get(peptide);
					if (log.motif!=null) {
						motifs.add(log.motif);
					}
				}
			}
			for (String motif : motifs) {
				System.out.println(motif);
			}
		}

		if (ANOVA_ANALYSIS) {
			for (int pep=0; pep<totalAdjustedPValues.length; pep++) {
				String peptide=peptides.get(pep);
				QuantitationLog log=primaryQuantLog.get(peptide);

				if (true||totalAdjustedPValues[pep]<0.05) { // NOTE REPORTS ALL, NOT JUST FDR CORRECTED

					float[][] data=log.getNormalizedData();
					double pValue=getANOVAPValue(data);

					System.out.println(log.peptideModSeq+"+"+log.charge+" rt:"+General.mean(log.rtInSecondsList.toArray())+" motif:"+log.motif+" localized:"+log.isSiteSpecific+" ("+log.protein+") p="+pValue+", FDR="+totalAdjustedPValues[pep]);

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
			System.out.print("Peptide\tProtein\tp-value\tFDR\tAKT\tLAKT\tMTOR\tMAPK");
			for (int i=0; i<6; i++) {
				System.out.print('\t');
				System.out.print(getSampleName(i+1));
			}
			System.out.println();

			ArrayList<String> flagged=new ArrayList<>();
			for (int pep=0; pep<totalAdjustedPValues.length; pep++) {
				if (totalAdjustedPValues[pep]<0.05) {
					String peptide=peptides.get(pep);
					QuantitationLog log=primaryQuantLog.get(peptide);

					float[][] data=log.getNormalizedData();
					double pValue=getANOVAPValue(data);
					System.out.print(log.peptideModSeq+"\t"+log.protein+"\t"+pValue+"\t"+totalAdjustedPValues[pep]);
					System.out.print("\t"+(log.motif!=null&&log.motif.matches("R.R..[ST].....")));
					System.out.print("\t"+(log.motif!=null&&log.motif.matches("..R..[ST].....")));
					System.out.print("\t"+(log.motif!=null&&log.motif.matches(".....[ST][FLW]....")));
					System.out.print("\t"+(log.motif!=null&&log.motif.matches(".....[ST]P....")));
					for (int i=0; i<data.length; i++) {
						System.out.print('\t');
						System.out.print(QuickMedian.median(data[i].clone()));
					}
					System.out.println();

					if (totalAdjustedPValues[pep]<pValue) {
						flagged.add(peptide);
					}
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
			HashMap<String, TDoubleArrayList> pvalueMap=new HashMap<>();
			HashMap<String, TFloatArrayList> rtMap=new HashMap<>();
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
			HashMap<String, double[]> pvalueMap=new HashMap<>();
			for (int pep=0; pep<totalAdjustedPValues.length; pep++) {
				String peptide=peptides.get(pep);
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

	private static Pair<HashMap<String, QuantitationLog>, HashMap<String, QuantitationLog>> getQuantData(String[] targets, File[] f) throws IOException, SQLException {
		HashMap<String, QuantitationLog> quantLog=new HashMap<>();
		HashMap<String, QuantitationLog> siteSpecificQuantLog=new HashMap<>();
		for (File file : f) {
			if (file.getName().endsWith(LibraryFile.ELIB)) {
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
							if (peptideSeq.indexOf(targets[i])>=0) {
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
		return new Pair<HashMap<String,QuantitationLog>, HashMap<String,QuantitationLog>>(quantLog, siteSpecificQuantLog);
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
				entry.motif=motif;
			}
		}
	}
	
	public static class QuantitationLog extends SimplePeptidePrecursor {
		boolean isSiteSpecific=false;
		final byte charge;
		String motif=null;
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
		public float[][] getNormalizedData() {
			float[][] data=getData();
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
