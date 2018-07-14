package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.apache.commons.math3.stat.inference.TestUtils;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SimplePeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.BenjaminiHochberg;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class P100ElibParser {
	private static final int numberOfSampleTypes=32;
	private static final int numberOfReplicates=3;
	private static String getSampleName(int i) {
		switch (i) {
		case 1: return "DMSO";
		case 2: return "Selumetinib";
		case 3: return "PD0325901";
		case 4: return "Everolimus";
		case 5: return "vemurafenib";
		case 6: return "TG101348";
		case 7: return "Tofacitinib";
		case 8: return "Pravastatin";
		case 9: return "flavopiridol";
		case 10: return "PD-0332991";
		case 11: return "Dinaciclib";
		case 12: return "RO4929097";
		case 13: return "BMS-906024";
		case 14: return "Verteporfin";
		case 15: return "SP600125";
		case 16: return "vorinostat";
		case 17: return "CC-401";
		case 18: return "SCH-900776";
		case 19: return "VX-970";
		case 20: return "losmapimod";
		case 21: return "PRI-724";
		case 22: return "dactolisib";
		case 23: return "afuresertib";
		case 24: return "BYL719";
		case 25: return "Pazopanib";
		case 26: return "Nilotinib";
		case 27: return "lenalidomide";
		case 28: return "AR-A014418";
		case 29: return "BMS-345541";
		case 30: return "IPI145";
		case 31: return "staurosporine";
		case 32: return "PS-1145";
		default: return i > 0 && i < 27 ? String.valueOf((char)(i + 64)) : null;
		}
	}
	public static TreeMap<String, Coordinate> sampleKey=new TreeMap<>();
	
	public static void main(String[] args) throws Exception {
		LibraryFile.OPEN_IN_PLACE=true;
		Logger.PRINT_TO_SCREEN=false;
		loadMap();

		String[] targets=null;//new String[] {"SFSKEVEER", "ILQEKLDQPVSAPPSPR", "HRGSEEDPLLSPVETWK", "RASGQAFELILSPR"};//KGSGDYMPMSPK;//targetPeptides;
		File[] f=new File("/Users/searleb/Documents/p100").listFiles();
		//File[] f=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/0.6.3_elibs").listFiles();
		
		//f=new File[] {new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/22jun2016_mcf7_phospho_1a.dia.thesaurus.elib")};
		
		Pair<TreeMap<String,QuantitationLog>, TreeMap<String,QuantitationLog>> quantLogPair=getQuantData(targets, f, new AminoAcidConstants());
		TreeMap<String, QuantitationLog> siteSpecificQuantLog=quantLogPair.y;
		
		int localized=0;
		ArrayList<String> siteSpecificPeptides=new ArrayList<>();
		TDoubleArrayList siteSpecificPValues=new TDoubleArrayList();
		for (String peptide : siteSpecificQuantLog.keySet()) {
			QuantitationLog log=siteSpecificQuantLog.get(peptide);
			//if (log.getNumMeasurements()<18) continue;
			if (!log.isAtLeastOneCaseFull()) continue;
			
			float[][] data=log.getNormalizedData();
			double pValue=getANOVAPValue(data);

			if (Double.isNaN(pValue)) continue;
			if (Double.isInfinite(pValue)) continue;
			if (pValue<0) continue;
			siteSpecificPeptides.add(peptide);
			siteSpecificPValues.add(pValue);
			if (log.bestLocalizationFDR<0.01) localized++;
		}
		
		System.out.println("Found "+siteSpecificPeptides.size()+"/"+siteSpecificQuantLog.size()+" peptides consistently ("+localized+" are site specific):");
		double[] siteSpecificAdjustedPValues=BenjaminiHochberg.calculateAdjustedPValues(siteSpecificPValues.toArray());

		System.out.println();
		System.out.print("\t\t\t\t\t\t\t");
		for (int i=0; i<numberOfSampleTypes; i++) {
			for (int j = 0; j < 3; j++) {
				System.out.print('\t');
				System.out.print(getSampleName(i+1)+(j+1));
			}
		}
		System.out.println();
		System.out.print("Peptide\tProtein\tlocalizationFDR\tquantFDR\tAKT\tLAKT\tMTOR\tMAPK");
		for (int i=0; i<numberOfSampleTypes; i++) {
			for (int j = 0; j < 3; j++) {
				System.out.print('\t');
				System.out.print(getSampleName(i+1));
			}
		}
		System.out.println();
		
		ArrayList<String> flagged=new ArrayList<>();
		for (int pep=0; pep<siteSpecificAdjustedPValues.length; pep++) {
			String peptide=siteSpecificPeptides.get(pep);
			QuantitationLog log=siteSpecificQuantLog.get(peptide);
			
			if (siteSpecificAdjustedPValues[pep]<0.05&&log.bestLocalizationFDR<0.01) {

				float[][] data=log.getNormalizedData();
				double pValue=getANOVAPValue(data);
				System.out.print(log.peptideModSeq+"\t"+log.toSitesString()+"\t"+log.bestLocalizationFDR+"\t"+siteSpecificAdjustedPValues[pep]);
				System.out.print("\t"+(log.doesMotifMatch("R.R..[ST].....")));
				System.out.print("\t"+(log.doesMotifMatch("..R..[ST].....")));
				System.out.print("\t"+(log.doesMotifMatch(".....[ST][FLW]....")));
				System.out.print("\t"+(log.doesMotifMatch(".....[ST]P....")));
				for (int i=0; i<data.length; i++) {
					for (int j = 0; j < data[i].length; j++) {
						System.out.print('\t');
						System.out.print(data[i][j]);
					}
				}
				System.out.println();

				if (siteSpecificAdjustedPValues[pep]<pValue) {
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

	private static double getANOVAPValue(float[][] data) {
		ArrayList<double[]> classes=new ArrayList<>();
		for (int i=0; i<data.length; i++) {
			classes.add(General.toDoubleArray(data[i]));
		}
		return TestUtils.oneWayAnovaPValue(classes);
	}
	
	
	private static void loadMap() {
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_A01_acq_01.mzML", new Coordinate(1,1));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_A02_acq_01.mzML", new Coordinate(2,1));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_A03_acq_01.mzML", new Coordinate(3,1));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_A04_acq_01.mzML", new Coordinate(1,9));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_A05_acq_01.mzML", new Coordinate(2,9));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_A06_acq_01.mzML", new Coordinate(3,9));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_A07_acq_01.mzML", new Coordinate(1,17));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_A08_acq_01.mzML", new Coordinate(2,17));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_A09_acq_01.mzML", new Coordinate(3,17));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_A10_acq_01.mzML", new Coordinate(1,25));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_A11_acq_01.mzML", new Coordinate(2,25));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_A12_acq_01.mzML", new Coordinate(3,25));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_B01_acq_01.mzML", new Coordinate(1,2));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_B02_acq_01.mzML", new Coordinate(2,2));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_B03_acq_01.mzML", new Coordinate(3,2));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_B04_acq_01.mzML", new Coordinate(1,10));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_B05_acq_01.mzML", new Coordinate(2,10));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_B06_acq_01.mzML", new Coordinate(3,10));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_B07_acq_01.mzML", new Coordinate(1,18));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_B08_acq_01.mzML", new Coordinate(2,18));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_B09_acq_01.mzML", new Coordinate(3,18));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_B10_acq_01.mzML", new Coordinate(1,26));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_B11_acq_01.mzML", new Coordinate(2,26));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_B12_acq_01.mzML", new Coordinate(3,26));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_C01_acq_01.mzML", new Coordinate(1,3));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_C02_acq_01.mzML", new Coordinate(2,3));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_C03_acq_01.mzML", new Coordinate(3,3));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_C04_acq_01.mzML", new Coordinate(1,11));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_C05_acq_01.mzML", new Coordinate(2,11));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_C06_acq_01.mzML", new Coordinate(3,11));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_C07_acq_01.mzML", new Coordinate(1,19));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_C08_acq_01.mzML", new Coordinate(2,19));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_C09_acq_01.mzML", new Coordinate(3,19));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_C10_acq_01.mzML", new Coordinate(1,27));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_C11_acq_01.mzML", new Coordinate(2,27));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_C12_acq_01.mzML", new Coordinate(3,27));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_D01_acq_01.mzML", new Coordinate(1,4));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_D02_acq_01.mzML", new Coordinate(2,4));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_D03_acq_01.mzML", new Coordinate(3,4));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_D04_acq_01.mzML", new Coordinate(1,12));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_D05_acq_01.mzML", new Coordinate(2,12));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_D06_acq_01.mzML", new Coordinate(3,12));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_D07_acq_01.mzML", new Coordinate(1,20));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_D08_acq_01.mzML", new Coordinate(2,20));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_D09_acq_01.mzML", new Coordinate(3,20));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_D10_acq_01.mzML", new Coordinate(1,28));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_D11_acq_01.mzML", new Coordinate(2,28));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_D12_acq_01.mzML", new Coordinate(3,28));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_E01_acq_01.mzML", new Coordinate(1,5));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_E02_acq_01.mzML", new Coordinate(2,5));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_E03_acq_01.mzML", new Coordinate(3,5));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_E04_acq_01.mzML", new Coordinate(1,13));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_E05_acq_01.mzML", new Coordinate(2,13));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_E06_acq_01.mzML", new Coordinate(3,13));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_E07_acq_01.mzML", new Coordinate(1,21));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_E08_acq_01.mzML", new Coordinate(2,21));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_E09_acq_01.mzML", new Coordinate(3,21));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_E10_acq_01.mzML", new Coordinate(1,29));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_E11_acq_01.mzML", new Coordinate(2,29));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_E12_acq_01.mzML", new Coordinate(3,29));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_F01_acq_01.mzML", new Coordinate(1,6));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_F02_acq_01.mzML", new Coordinate(2,6));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_F03_acq_01.mzML", new Coordinate(3,6));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_F04_acq_01.mzML", new Coordinate(1,14));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_F05_acq_01.mzML", new Coordinate(2,14));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_F06_acq_01.mzML", new Coordinate(3,14));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_F07_acq_01.mzML", new Coordinate(1,22));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_F08_acq_01.mzML", new Coordinate(2,22));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_F09_acq_01.mzML", new Coordinate(3,22));
		sampleKey.put("CC20160706_P100_Plate34_PC3_T3_P-0034_F10_acq_01.mzML", new Coordinate(1,30));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_F11_acq_02.mzML", new Coordinate(2,30));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_F12_acq_02.mzML", new Coordinate(3,30));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_G01_acq_01.mzML", new Coordinate(1,7));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_G02_acq_01.mzML", new Coordinate(2,7));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_G03_acq_01.mzML", new Coordinate(3,7));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_G04_acq_01.mzML", new Coordinate(1,15));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_G05_acq_01.mzML", new Coordinate(2,15));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_G06_acq_01.mzML", new Coordinate(3,15));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_G07_acq_01.mzML", new Coordinate(1,23));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_G08_acq_01.mzML", new Coordinate(2,23));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_G09_acq_01.mzML", new Coordinate(3,23));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_G10_acq_01.mzML", new Coordinate(1,31));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_G11_acq_02.mzML", new Coordinate(2,31));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_G12_acq_01.mzML", new Coordinate(3,31));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_H01_acq_01.mzML", new Coordinate(1,8));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_H02_acq_01.mzML", new Coordinate(2,8));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_H03_acq_01.mzML", new Coordinate(3,8));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_H04_acq_01.mzML", new Coordinate(1,16));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_H05_acq_01.mzML", new Coordinate(2,16));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_H06_acq_01.mzML", new Coordinate(3,16));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_H07_acq_01.mzML", new Coordinate(1,24));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_H08_acq_01.mzML", new Coordinate(2,24));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_H09_acq_01.mzML", new Coordinate(3,24));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_H10_acq_01.mzML", new Coordinate(1,32));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_H11_acq_01.mzML", new Coordinate(2,32));
		sampleKey.put("CN20160706_P100_Plate34_PC3_T3_P-0034_H12_acq_01.mzML", new Coordinate(3,32));
	}


	private static Pair<TreeMap<String, QuantitationLog>, TreeMap<String, QuantitationLog>> getQuantData(String[] targets, File[] f, AminoAcidConstants aaConstants) throws IOException, SQLException {
		TreeMap<String, QuantitationLog> quantLog=new TreeMap<>();
		TreeMap<String, QuantitationLog> siteSpecificQuantLog=new TreeMap<>();
		for (File file : f) {
			if (file.getName().endsWith(".thesaurus.elib")) {
				System.out.println("Parsing "+file.getName()+"...");
				LibraryFile library=new LibraryFile();
				library.openFile(file);
				

				Connection c=library.getConnection();
				Statement s=c.createStatement();
				ResultSet rs = s.executeQuery("select pep.PrecursorCharge, pep.PeptideModSeq, pep.PeptideSeq, pep.SourceFile, max(pep.LocalizedIntensity), max(pep.TotalIntensity), pep.localizationFDR, pep.RTInSecondsCenter,pep.localizationScore,"+
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
					double bestLocalizationFDR=rs.getDouble(7);
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
							log=new QuantitationLog(proteinToken, peptideModSeq, precursorCharge, aaConstants);
							quantLog.put(peptideModSeq, log);
							siteLog=new QuantitationLog(proteinToken, peptideModSeq, precursorCharge, aaConstants);
							siteSpecificQuantLog.put(peptideModSeq, siteLog);
						}
						log.addIntensity(coord, totalIntensity, rtInSeconds, localizationScore, bestLocalizationFDR);
						siteLog.addIntensity(coord, localizedIntensity, rtInSeconds, localizationScore, bestLocalizationFDR);
						
					}
				}
				rs.close();
				s.close();
				c.close();
			}
		}
		return new Pair<TreeMap<String,QuantitationLog>, TreeMap<String,QuantitationLog>>(quantLog, siteSpecificQuantLog);
	}
	

	public static class QuantitationLog extends SimplePeptidePrecursor {
		double bestLocalizationFDR=1.0;
		final byte charge;
		TreeMap<String, TIntObjectHashMap<String>> motifMap=new TreeMap<String, TIntObjectHashMap<String>>();
		final String protein;
		final String peptideModSeq;
		final TObjectFloatHashMap<Coordinate> intensities=new TObjectFloatHashMap<>();
		final TFloatArrayList rtInSecondsList=new TFloatArrayList();
		final TFloatArrayList localizationScores=new TFloatArrayList();
		
		public QuantitationLog(String protein, String peptideModSeq, byte charge, AminoAcidConstants aaConstants) {
			super(peptideModSeq, charge, aaConstants);
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
		
		public void addIntensity(Coordinate c, float intensity, float rtInSeconds, float localizationScore, double bestLocalizationFDR) {
			intensities.adjustOrPutValue(c, intensity, intensity);
			if (this.bestLocalizationFDR>bestLocalizationFDR) {
				this.bestLocalizationFDR=bestLocalizationFDR;
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
