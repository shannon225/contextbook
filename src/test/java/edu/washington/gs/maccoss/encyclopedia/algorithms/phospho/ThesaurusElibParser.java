package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.stat.inference.TestUtils;

import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.BenjaminiHochberg;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class ThesaurusElibParser {
	public static HashMap<String, Coordinate> sampleKey=new HashMap<>();
	
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
	
	public static void main(String[] args) throws Exception {
		StripeFile.OPEN_IN_PLACE=true;
		LibraryFile.OPEN_IN_PLACE=true;
		Logger.PRINT_TO_SCREEN=false;
		loadMap();
		
		File[] f=new File("/Users/searleb/Documents/school/localization_manuscript/elibs/mcf7").listFiles();
		
		HashMap<String, QuantitationLog> quantLog=new HashMap<>();
		for (File file : f) {
			if (file.getName().endsWith(LibraryFile.ELIB)) {
				System.out.println("Parsing "+file.getName()+"...");
				LibraryFile library=new LibraryFile();
				library.openFile(file);
				

				Connection c=library.getConnection();
				Statement s=c.createStatement();
				ResultSet rs=s.executeQuery("select pep.PrecursorCharge, pep.PeptideModSeq, pep.PeptideSeq, pep.SourceFile, pep.TotalIntensity, pep.IsSiteSpecific, pro.ProteinAccessions from peptidelocalizations pep, proteins pro where pep.PeptideSeq = pro.PeptideSeq");

				while (rs.next()) {
					//byte precursorCharge=(byte)rs.getInt(1);
					String peptideModSeq=rs.getString(2);
					String peptideSeq=rs.getString(3);
					String sourceFile=rs.getString(4);
					float totalIntensity=rs.getFloat(5);
					boolean isSiteSpecific=rs.getBoolean(6);
					String proteinToken=rs.getString(7);
					//HashSet<String> accessions=PSMData.stringToAccessions(proteinToken);
					boolean keeper=false;
					for (int i=0; i<targetPeptides.length; i++) {
						if (peptideSeq.indexOf(targetPeptides[i])>=0) {
							keeper=true;
							break;
						}
					}
					
					if (keeper) {
						Coordinate coord=sampleKey.get(sourceFile);
						if (coord==null) {
							System.out.println("FAILED TO FIND SAMPLE: "+sourceFile);
							System.exit(1);
						}

						QuantitationLog log=quantLog.get(peptideModSeq);
						if (log==null) {
							log=new QuantitationLog(proteinToken, peptideModSeq);
							quantLog.put(peptideModSeq, log);
						}
						log.addIntensity(coord, totalIntensity, isSiteSpecific);
					}
				}
				rs.close();
				s.close();
				c.close();
			}
		}
		
		ArrayList<String> peptides=new ArrayList<>();
		TDoubleArrayList pValues=new TDoubleArrayList();
		for (String peptide : quantLog.keySet()) {
			QuantitationLog log=quantLog.get(peptide);
			if (log.getNumMeasurements()<16) continue;
			
			float[][] data=log.getNormalizedData();
			double pValue=getPValue(data);
			
			if (Double.isNaN(pValue)) continue;
			if (Double.isInfinite(pValue)) continue;
			peptides.add(peptide);
			pValues.add(pValue);
		}
		
		double[] adjustedPValues=BenjaminiHochberg.calculateAdjustedPValues(pValues.toArray());
		for (int pep=0; pep<adjustedPValues.length; pep++) {
			if (adjustedPValues[pep]<0.05) {
				String peptide=peptides.get(pep);
				QuantitationLog log=quantLog.get(peptide);

				float[][] data=log.getNormalizedData();
				double pValue=getPValue(data);

				System.out.println(log.peptide+" localized:"+log.isSiteSpecific+" ("+log.protein+") p="+pValue+", FDR="+adjustedPValues[pep]);

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

		System.out.println();
		System.out.print("Peptide\tProtein");
		for (int i=0; i<6; i++) {
			System.out.print('\t');
			System.out.print(getSampleName(i+1));
		}
		System.out.println();
		
		for (int pep=0; pep<adjustedPValues.length; pep++) {
			if (adjustedPValues[pep]<0.05) {
				String peptide=peptides.get(pep);
				QuantitationLog log=quantLog.get(peptide);

				float[][] data=log.getNormalizedData();
				System.out.print(log.peptide+"\t"+log.protein);
				for (int i=0; i<data.length; i++) {
					System.out.print('\t');
					System.out.print(QuickMedian.median(data[i].clone()));
				}
				System.out.println();
			}
		}
	}
	
	private static double getPValue(float[][] data) {
		ArrayList<double[]> classes=new ArrayList<>();
		for (int i=0; i<data.length; i++) {
			classes.add(General.toDoubleArray(data[i]));
		}
		return TestUtils.oneWayAnovaPValue(classes);
	}
	
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
	
	public static class QuantitationLog {
		boolean isSiteSpecific=false;
		final String protein;
		final String peptide;
		final TObjectFloatHashMap<Coordinate> intensities=new TObjectFloatHashMap<>();
		
		public QuantitationLog(String protein, String peptide) {
			this.protein=protein;
			this.peptide=peptide;
		}
		
		public int getNumMeasurements() {
			return intensities.size();
		}
		
		public void addIntensity(Coordinate c, float intensity, boolean isSiteSpecific) {
			intensities.adjustOrPutValue(c, intensity, intensity);
			if (isSiteSpecific) {
				this.isSiteSpecific=true;
			}
		}
		public float[][] getData() {
			Coordinate[] coords=intensities.keys(new Coordinate[intensities.size()]);
			int maxReplicate=0;
			int maxSample=0;
			for (Coordinate c : coords) {
				if (c.replicate>maxReplicate) maxReplicate=c.replicate;
				if (c.sample>maxSample) maxSample=c.sample;
			}
			
			float[][] results=new float[maxSample][];
			for (int i=0; i<results.length; i++) {
				results[i]=new float[maxReplicate];
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
