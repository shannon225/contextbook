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
import java.util.TreeMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeSet;

import org.apache.commons.math3.stat.inference.TestUtils;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideTrie;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SimplePeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.StringUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.BenjaminiHochberg;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class ThesaurusElibParser {
	public static TreeMap<String, Coordinate> sampleKey=new TreeMap<>();
	
	private static final int numberOfSampleTypes=6;
	private static final int numberOfReplicates=6;
	private static String getSampleName(int i) {
		switch (i) {
		case 1: return "Cont";
		case 2: return "Ins";
		case 3: return "IGF1";
		case 4: return "MK-Cont";
		case 5: return "MK-Ins";
		case 6: return "MK-IGF1";
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
	
	private static final byte INS_VS_CONTROL=1;
	private static final byte IGF_VS_CONTROL=2;
	private static final byte INS_VS_IGF=3;
	private static final byte INS_IGF_VS_MK2206=4;
	private static final double getPairedTTest(float[][] data, byte testType) {
		TDoubleArrayList x=new TDoubleArrayList();
		TDoubleArrayList y=new TDoubleArrayList();
		
		switch (testType) {
		case INS_VS_CONTROL:
			x.addAll(General.toDoubleArray(data[0]));
			x.addAll(General.toDoubleArray(data[3]));
			y.addAll(General.toDoubleArray(data[1]));
			y.addAll(General.toDoubleArray(data[4]));
			break;
		case IGF_VS_CONTROL:
			x.addAll(General.toDoubleArray(data[0]));
			x.addAll(General.toDoubleArray(data[3]));
			y.addAll(General.toDoubleArray(data[2]));
			y.addAll(General.toDoubleArray(data[5]));
			break;
		case INS_VS_IGF:
			x.addAll(General.toDoubleArray(data[1]));
			x.addAll(General.toDoubleArray(data[4]));
			y.addAll(General.toDoubleArray(data[2]));
			y.addAll(General.toDoubleArray(data[5]));
			break;
		case INS_IGF_VS_MK2206:
			x.addAll(General.toDoubleArray(data[1]));
			x.addAll(General.toDoubleArray(data[2]));
			y.addAll(General.toDoubleArray(data[4]));
			y.addAll(General.toDoubleArray(data[5]));
			break;
		default:
			break;
		}
		
		return TestUtils.pairedTTest(x.toArray(), y.toArray());
	}
	
	
	public static String[] targetProteins=new String[] { "O43521", "O43524", "O60343", "O60825", "O75581", "P02545", "P04049", "P04637", "P06239", "P06730", "P07948", "P08069", "P10415", "P11274",
			"P12931", "P13807", "P22681", "P23443", "P27361", "P29474", "P31749", "P31751", "P35568", "P42345", "P43403", "P45983", "P45984", "P49023", "P49815", "P49840", "P49841", "P51812",
			"P53396", "P54646", "P62136", "P62753", "P98177", "Q00987", "Q02750", "Q03135", "Q05397", "Q12778", "Q13131", "Q13164", "Q13322", "Q13480", "Q13541", "Q15418", "Q6R327", "Q96B36",
			"Q96BR1", "Q9UQC2" };
	
	public static String[] targetPeptides=new String[] { "AALQTAPESADDSPSQLSKWPGSPTSR", "AASMDNNSK", "AASMDSSSK", "AEERPTFDYLQSVLDDFYTATEGQYQQQP", "AENGLLMTPCYTANFVAPEVLK", "AENGLLMTPCYTANFVAPEVLKR",
			"AFPEHFTYEPNEADAAQGYRYPRPASVPPSPSLSR", "AHTFSHPPSSTK", "AHTFSHPPSSTKR", "AISETEENSDELSGER", "AISETEENSDELSGERQR", "ALPNNTSSSPQPK", "ARTSSFAEPGGGGGGGGGGPGGSASGPGGTGGGK", "AVSMDNSNK",
			"AVSMDNSNKYTK", "AYSFCGTVEYMAPEVVNRR", "CSSVTGVQR", "DGATMKTFCGTPEYLAPEVLEDNDYGR", "DIYETDYYR", "DIYETDYYRK", "DIYETDYYRKGGK", "EGIAISDTTTTFCGTPEYLAPEVIR",
			"EGIKDGATMKTFCGTPEYLAPEVLEDNDYGR", "EGISDGATMKTFCGTPEYLAPEVLEDNDYGR", "EKAEERPTFDYLQSVLDDFYTATEGQYQQQP", "ELELMFGCQVEGDAAETPPRPR", "FIGSPRTPVSPVK", "FLMECRNSPVTK", "FLMECRNSPVTKTPPR",
			"FTRQTPVDSPDDSTLSESANQVFLGFTYVAPSVLESVK", "FTRQTPVDSPDDSTLSESANQVFLGFTYVAPSVLESVKEK", "GDKQVEYLDLDLDSGK", "GFSFVATGLMEDDGKPR", "GFSFVATGLMEDDGKPRAPQAPLHSVVQQLHGK",
			"GGHHRPDSSTLHTDDGYMPMSPGVAPVPSGR", "GHGQPGADAEKPFYVNVEFHHER", "GKYGQFSGLNPGGRPITPPR", "GKYGQFSGLNPGGRPITPPRNSAK", "GLCTSPAEHQYFMTEYVATR", "GRLGSVDSFER", "GSGDYMPMSPK",
			"GTYFPAILNPPPSPATER", "GYTISDSAPSR", "GYTISDSAPSRR", "HSSETFSSTPSATR", "HVSISYDIPPTPGNTYQIPR", "IADFGLSNMMSDGEFLRTSCGSPNYAAPEVISGR", "IADPEHDHTGFLTEYVATR", "IGDFGMTRDIYETDYYRK",
			"IQAAASTPTNATAASDANTGDR", "IQAAASTPTNATAASDANTGDRGQTNNAASASASNST", "IRTLTEPSVDFNHSDDFTPISTVQK", "ITSPDKPRPMPMDTSVYESPYSDPEELK", "ITSPDKPRPMPMDTSVYESPYSDPEELKDK", "IVIGYQSHADTATKSGSTTK",
			"IVIGYQSHADTATKSGSTTKNR", "KAYSFCGTVEYMAPEVVNRR", "KFLMECRNSPVTK", "KGSGDYMPMSPK", "KQEEEEMDFRSGSPSDNSGAEEMEVSLAKPK", "KTGTTVPESIHSFIGDGLVKPEALNK", "KTGTTVPESIHSFIGDGLVKPEALNKK",
			"LCDFGVSGQLIDSMANSFVGTR", "LGSVDSFER", "LIEDNEYTAR", "LIEDNEYTAREGAK", "LIEDNEYTAREGAKFPIK", "LIEDNEYTARQGAK", "LMFKTEGPDSD", "LNTSDFQK", "LNTSDFQKLK", "LPPGEQCEGEEDTEYMTPSSRPLRPLDTSQSSR",
			"LSSLRASTSKSESSQK", "MEEPQSDPSVEPPLSQETFSDLWK", "MNILGSQSPLHPSTLSTVIHR", "NSPVTKTPPR", "NSPVTKTPPRDLPTIPGVTSPSSDEPPMEASQSHLR", "NYSVGSRPLKPLSPLR", "QLRAENGLLMTPCYTANFVAPEVLK",
			"QLRAENGLLMTPCYTANFVAPEVLKR", "QTPVDSPDDSTLSESANQVFLGFTYVAPSVLESVKEK", "QVEYLDLDLDSGK", "RAHTFSHPPSSTK", "RAHTFSHPPSSTKR", "RAISETEENSDELSGER", "RAISETEENSDELSGERQR", "RALPNNTSSSPQPK",
			"RAVSMDNSNKYTK", "RFIGSPRTPVSPVK", "RGGHHRPDSSTLHTDDGYMPMSPGVAPVPSGR", "RHSSETFSSTPSATR", "RLSSLRASTSK", "RPHFPQFSYSASGTA", "RRAISETEENSDELSGER", "RRHSSETFSSTPSATR",
			"RVVLGDGVQLPPGDYSTTPGGTLFSTTPGGTR", "SDSTNSEDNYVPMNPGSSTLLAMER", "SGAQASSTPLSPTR", "SGAQASSTPLSPTRITR", "SGSPSDNSGAEEMEVSLAKPK", "SGTATPQR", "SHSESASPSALSSSPNNLSPTGWSQPK",
			"SHSESASPSALSSSPNNLSPTGWSQPKTPVPAQR", "SIDDEITEAKSGTATPQR", "SLPVSVPVWGFK", "SLPVSVPVWGFKEK", "SPGEYVNIEFGSDQSGYLSGPVAFHSSPSVR", "SPLFIFMR", "SRCSSVTGVQR", "SRTESITATSPASMVGGKPGSFR",
			"SSGSGSSVADERVDYVVVDQQK", "SSSFPYTTK", "SSSSPELQTLQDILGDPGDK", "SSSSPELQTLQDILGDPGDKADVGR", "STSLNERPK", "STSLNERPKR", "STSTPNVHMVSTTLPVDSR", "TACTNFMMTPYVVTR", "TAGTSFMMTPYVVTR",
			"TASFSESR", "TASFSESRADEVAPAK", "TASFSESRADEVAPAKK", "TDSYSAGQSVEILDGVELGEPAHK", "TDSYSAGQSVEILDGVELGEPAHKK", "TESITATSPASMVGGKPGSFR", "TFCGTPEYLAPEVLEDNDYGR", "TGTTVPESIHSFIGDGLVKPEALNK",
			"TGTTVPESIHSFIGDGLVKPEALNKK", "THAVSVSETDDYAEIIDEEDTYTMPSTR", "THFPQFSYSASIRE", "TLTEPSVDFNHSDDFTPISTVQK", "TPGRPLSSYGMDSRPPMAIFELLDYIVNEPPPK", "TPKDSPGIPPSAGAHQLFR",
			"TPPRDLPTIPGVTSPSSDEPPMEASQSHLR", "TPVSPVKFSPGDFWGR", "TQSFSLQER", "TRTDSYSAGQSVEILDGVELGEPAHK", "TRTDSYSAGQSVEILDGVELGEPAHKK", "TSPLQTPAAPGAAAGPALSPVPPVVHLTLR",
			"TSSFAEPGGGGGGGGGGPGGSASGPGGTGGGK", "TSSNASTISGRLSPIMTEQDDLGEGDVHSMVYPPSAAK", "TSSVSNPQDSVGSPCSRVGEEEHVYSFPNK", "TTSFAESCKPVQQPSAFGSMK", "VDYVVVDQQK", "VGEEEHVYSFPNK", "VGEEEHVYSFPNKQK",
			"VIEDNEYTAR", "VIEDNEYTAREGAK", "VVLGDGVQLPPGDYSTTPGGTLFSTTPGGTR", "VYENVTGLVK", "WPGSPTSR", "YFDDEFTAQSITITPPDR", "YFDDEFTAQSITITPPDRYDSLGLLELDQR", "YGQFSGLNPGGRPITPPR",
			"YGQFSGLNPGGRPITPPRNSAK", "YMEDSTYYK", "YMEDSTYYKASK", "YPRPASVPPSPSLSR", "YVDSEGHLYTVPIR", "YVDSEGHLYTVPIREQGNIYKPNNK" };
	
	public static String[] KGSGDYMPMSPK=new String[] {"KGSGDYMPMSPK"};

	public static void loadMap() {
		Arrays.sort(targetProteins);
		
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
		sampleKey.put("22jun2016_mcf7_phospho_5a.mzML", new Coordinate(5, 1));
		sampleKey.put("22jun2016_mcf7_phospho_5b.mzML", new Coordinate(5, 2));
		sampleKey.put("22jun2016_mcf7_phospho_5c.mzML", new Coordinate(5, 3));
		sampleKey.put("22jun2016_mcf7_phospho_5d.mzML", new Coordinate(5, 4));
		sampleKey.put("22jun2016_mcf7_phospho_5e.mzML", new Coordinate(5, 5));
		sampleKey.put("22jun2016_mcf7_phospho_5f.mzML", new Coordinate(5, 6));
		sampleKey.put("22jun2016_mcf7_phospho_6a.mzML", new Coordinate(6, 1));
		sampleKey.put("22jun2016_mcf7_phospho_6b.mzML", new Coordinate(6, 2));
		sampleKey.put("22jun2016_mcf7_phospho_6c.mzML", new Coordinate(6, 3));
		sampleKey.put("22jun2016_mcf7_phospho_6d.mzML", new Coordinate(6, 4));
		sampleKey.put("22jun2016_mcf7_phospho_6e.mzML", new Coordinate(6, 5));
		sampleKey.put("22jun2016_mcf7_phospho_6f.mzML", new Coordinate(6, 6));
	}
	public static final boolean TOTAL_ANALYSIS=true;
	
	public static final boolean MOTIF_ANALYSIS=false;
	public static final boolean ANOVA_ANALYSIS=true;
	public static final boolean HEATMAP_ANALYSIS=false;
	public static final boolean MULTIPLE_FORM_ANALYSIS=false;
	public static final boolean SITE_SPECIFIC_VS_TOTAL_ANALYSIS=false;
	
	public static void main(String[] args) throws Exception {
		StripeFile.OPEN_IN_PLACE=true;
		LibraryFile.OPEN_IN_PLACE=true;
		Logger.PRINT_TO_SCREEN=false;
		loadMap();

		PeptideModification mod=PeptideModification.phosphorylation;
		String[] targets=null;//new String[] {"SFSKEVEER", "ILQEKLDQPVSAPPSPR", "HRGSEEDPLLSPVETWK", "RASGQAFELILSPR"};//KGSGDYMPMSPK;//targetPeptides;
		File[] f=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs").listFiles();
		//f=new File[] {new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/22jun2016_mcf7_phospho_1a.dia.thesaurus.elib")};
		
		Pair<TreeMap<String,QuantitationLog>, TreeMap<String,QuantitationLog>> quantLogPair=getQuantData(targets, f);
		TreeMap<String, QuantitationLog> totalQuantLog=quantLogPair.x;
		TreeMap<String, QuantitationLog> siteSpecificQuantLog=quantLogPair.y;
		
		TreeMap<String, QuantitationLog> primaryQuantLog=TOTAL_ANALYSIS?totalQuantLog:siteSpecificQuantLog;
		
		PeptideMotifTrie motifTrie=new PeptideMotifTrie(primaryQuantLog.values(), mod);

		System.out.println("Reading FASTA...");
		ArrayList<FastaEntryInterface> fasta=FastaReader.readFasta(new File("/Users/searleb/Downloads/uniprot.HUMAN.fasta"));
		//ArrayList<FastaEntryInterface> fasta=FastaReader.readFasta(new File("/Users/searleb/Documents/school/projects/pecandata/UP000005640_9606.fasta"));
		motifTrie.addFasta(fasta);
		
		ArrayList<String> siteSpecificPeptides=new ArrayList<>();
		TDoubleArrayList siteSpecificPValues=new TDoubleArrayList();
		int localized=0;
		for (String peptide : siteSpecificQuantLog.keySet()) {
			QuantitationLog log=siteSpecificQuantLog.get(peptide);
			//if (log.getNumMeasurements()<18) continue;
			if (!log.isAtLeastOneCaseFull()) continue;
			
			float[][] data=log.getNormalizedData();
			double pValue;
			
			pValue=getANOVAPValue(data);
			//pValue=General.max(new double[] {getPairedTTest(data, INS_VS_CONTROL), getPairedTTest(data, IGF_VS_CONTROL), getPairedTTest(data, INS_IGF_VS_MK2206)}); // MK2206 changes!
			//pValue=General.max(new double[] {General.min(new double[] {getPairedTTest(data, INS_VS_CONTROL), getPairedTTest(data, IGF_VS_CONTROL)})}); // INS and IGF
			
			if (Double.isNaN(pValue)) continue;
			if (Double.isInfinite(pValue)) continue;
			if (pValue<0) continue;
			siteSpecificPeptides.add(peptide);
			siteSpecificPValues.add(pValue);
			if (log.isSiteSpecific) localized++;
		}
		System.out.println("Found "+siteSpecificPeptides.size()+" peptides ("+localized+" are site specific):");
		double[] siteSpecificAdjustedPValues=BenjaminiHochberg.calculateAdjustedPValues(siteSpecificPValues.toArray());

		ArrayList<String> totalPeptides=new ArrayList<>();
		TDoubleArrayList totalPValues=new TDoubleArrayList();
		for (String peptide : totalQuantLog.keySet()) {
			QuantitationLog log=totalQuantLog.get(peptide);
			//if (log.getNumMeasurements()<18) continue;
			if (!log.isAtLeastOneCaseFull()) continue;
			
			float[][] data=log.getNormalizedData();
			double pValue;
			
			pValue=getANOVAPValue(data);
			//pValue=General.max(new double[] {getPairedTTest(data, INS_VS_CONTROL), getPairedTTest(data, IGF_VS_CONTROL), getPairedTTest(data, INS_IGF_VS_MK2206)}); // MK2206 changes!
			//pValue=General.max(new double[] {General.min(new double[] {getPairedTTest(data, INS_VS_CONTROL), getPairedTTest(data, IGF_VS_CONTROL)})}); // INS and IGF
			
			if (Double.isNaN(pValue)) continue;
			if (Double.isInfinite(pValue)) continue;
			if (pValue<0) continue;
			totalPeptides.add(peptide);
			totalPValues.add(pValue);
		}
		double[] totalAdjustedPValues=BenjaminiHochberg.calculateAdjustedPValues(totalPValues.toArray());
		
		if (MOTIF_ANALYSIS) {
			TreeSet<String> motifs=new TreeSet<>();

			for (int pep=0; pep<totalAdjustedPValues.length; pep++) {
				if (totalAdjustedPValues[pep]<0.05) {
					String peptide=totalPeptides.get(pep);
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
			int length=primaryQuantLog==siteSpecificQuantLog?siteSpecificPValues.size():totalPValues.size();
			for (int pep=0; pep<length; pep++) {
				// if
				// (log.motif!=null&&log.motif.matches("R.R..[ST].....")&&adjustedPValues[pep]<0.05)
				// {
				
				String peptide;
				QuantitationLog log;
				double fdr;
				double pValue;
				
				if (primaryQuantLog==siteSpecificQuantLog) {
					peptide=siteSpecificPeptides.get(pep);
					log=siteSpecificQuantLog.get(peptide);
					pValue=siteSpecificPValues.get(pep);
					fdr=siteSpecificAdjustedPValues[pep];
				} else {
					peptide=totalPeptides.get(pep);
					log=totalQuantLog.get(peptide);
					pValue=totalPValues.get(pep);
					fdr=totalAdjustedPValues[pep];
				}
				if (true||fdr<0.05) {

					float[][] data=log.getNormalizedData();

					System.out.println(log.peptideModSeq+"+"+log.charge+" rt:"+General.mean(log.rtInSecondsList.toArray())+" localized:"+log.isSiteSpecific+" p="+pValue+", FDR="+fdr+" ("+log.toSitesString()+")");
					
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
							System.out.print(data[samp][rep]);
						}
						System.out.println();
					}
					System.out.println();
				}
			}
		}

		if (HEATMAP_ANALYSIS) {
			System.out.println();
			System.out.print("Peptide\tProtein\tisSiteSpecific\tp-value\tFDR\tAKT\tLAKT\tMTOR\tMAPK");
			for (int i=0; i<6; i++) {
				System.out.print('\t');
				System.out.print(getSampleName(i+1));
			}
			System.out.println();

			ArrayList<String> flagged=new ArrayList<>();
			for (int pep=0; pep<totalAdjustedPValues.length; pep++) {
				if (totalAdjustedPValues[pep]<0.05) {
					String peptide=totalPeptides.get(pep);
					QuantitationLog log=primaryQuantLog.get(peptide);

					float[][] data=log.getNormalizedData();
					double pValue=getANOVAPValue(data);
					System.out.print(log.peptideModSeq+"\t"+log.toSitesString()+"\t"+log.isSiteSpecific+"\t"+pValue+"\t"+totalAdjustedPValues[pep]);
					System.out.print("\t"+(log.doesMotifMatch("R.R..[ST].....")));
					System.out.print("\t"+(log.doesMotifMatch("..R..[ST].....")));
					System.out.print("\t"+(log.doesMotifMatch(".....[ST][FLW]....")));
					System.out.print("\t"+(log.doesMotifMatch(".....[ST]P....")));
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
			TreeMap<String, TDoubleArrayList> pvalueMap=new TreeMap<>();
			TreeMap<String, TFloatArrayList> rtMap=new TreeMap<>();
			for (int pep=0; pep<totalAdjustedPValues.length; pep++) {
				String peptide=totalPeptides.get(pep);
				QuantitationLog log=totalQuantLog.get(peptide);
				String key=peptide.replace("[+79.966331]", "");
				
				TDoubleArrayList list=pvalueMap.get(key);
				TFloatArrayList rtList=rtMap.get(key);
				if (list==null) {
					list=new TDoubleArrayList();
					pvalueMap.put(key, list);
					rtList=new TFloatArrayList();
					rtMap.put(key, rtList);
				}
				float control = Math.max(1, QuickMedian.median(log.getData()[0].clone()));
				float insulin = Math.max(1, QuickMedian.median(log.getData()[1].clone()));
				list.add(Log.log2(insulin/control));
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
				String peptide=totalPeptides.get(pep);
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

	private static Pair<TreeMap<String, QuantitationLog>, TreeMap<String, QuantitationLog>> getQuantData(String[] targets, File[] f) throws IOException, SQLException {
		TreeMap<String, QuantitationLog> quantLog=new TreeMap<>();
		TreeMap<String, QuantitationLog> siteSpecificQuantLog=new TreeMap<>();
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
		return new Pair<TreeMap<String,QuantitationLog>, TreeMap<String,QuantitationLog>>(quantLog, siteSpecificQuantLog);
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
