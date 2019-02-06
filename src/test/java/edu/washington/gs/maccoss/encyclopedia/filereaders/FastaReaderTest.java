package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;

import org.apache.commons.math3.util.CombinatoricsUtils;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.ExtendedFastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filewriters.FastaWriter;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

public class FastaReaderTest extends TestCase {
	/*public static void main(String[] args) throws IOException, FileNotFoundException {
		File f=new File("/Users/searleb/Downloads/hg38_6FT.fasta");
		File out=new File("/Users/searleb/Downloads/hg38_coding.fasta");
		BufferedReader in=new BufferedReader(new FileReader(f));
		FastaWriter writer=new FastaWriter(out);

		boolean inAnnotation=false;
		StringBuilder annotation=new StringBuilder();
		StringBuilder sequence=new StringBuilder();
		char[] buffer=new char[1024*1024];
		int length;
		int counter=0;
		int index=0;
		int startIndex=0;
		while ((length=in.read(buffer))>=0) {
			counter++;
			if (counter%100==0) {
				System.out.println(" "+counter+" MB");
			} else if (counter%10==0) {
				System.out.print(". ");
			} else {
				System.out.print('.');
			}
			
			for (int i=0; i<length; i++) {
				index++;
				if (buffer[i]=='>') {
					inAnnotation=true;
					annotation.setLength(0);
				} else if (buffer[i]=='\n') {
					if (sequence.length()>=8) {
						writer.write(new FastaEntry(f.getName(), annotation.toString()+"."+(startIndex+1)+"."+(index-1), sequence.toString()));
					}
					sequence.setLength(0);
					
					inAnnotation=false;
					index=0;
					startIndex=index;
				} else {
					if (inAnnotation) {
						annotation.append(buffer[i]);
					} else {
						if (buffer[i]=='*'||buffer[i]=='X') {
							if (sequence.length()>=8) {
								writer.write(new FastaEntry(f.getName(), annotation.toString()+"."+(startIndex+1)+"."+(index-1), sequence.toString()));
							}
							sequence.setLength(0);
							startIndex=index;
						} else {
							sequence.append(buffer[i]);
						}
					}
				}
			}
		}
		in.close();
		writer.close();
	}*/
	
	
	/**
	 * for pecan peptide listing
	 * @param args
	 * @throws Exception
	 */
	public static void main2(String[] args) throws Exception {
		PecanSearchParameters parameters=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false, true, false);
		//File f=new File("/Users/searleb/Documents/projects/phosphopedia/sp_iso_HUMAN_4.9.2015_UP000005640.fasta");
		File f=new File("/Users/searleb/Documents/chromatogram_library_manuscript/real_pecan/cerevisiae_orf_trans_all.fasta");
		ArrayList<FastaEntryInterface> targetProteins=FastaReader.readFasta(f, parameters);
		
		//PrintWriter writer=new PrintWriter("/Users/searleb/Documents/chromatogram_library_manuscript/sp_iso_HUMAN_4.9.2015_UP000005640.peptides.txt");
		PrintWriter writer=new PrintWriter("/Users/searleb/Documents/chromatogram_library_manuscript/real_pecan/cerevisiae_orf_trans_all.peptides.txt");
		for (FastaEntryInterface entry : targetProteins) {
			ArrayList<FastaPeptideEntry> peptides=parameters.getEnzyme().digestProtein(entry, parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages(), parameters.getAAConstants(), false);
			for (FastaPeptideEntry peptide : peptides) {
				writer.println(entry.getAccession()+"\t"+PeptideUtils.getPeptideSeq(peptide.getSequence()));
			}
		}
		writer.flush();
		writer.close();
	}
	
	public static void mainx(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		File f=new File("/Users/bsearle/Documents/prosit/hela/uniprot-9606.fasta");
		
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		
		for (FastaEntryInterface entry : entries) {
			System.out.println(entry.getSequence().length());
		}
	}
	
	public static void main(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		//File f=new File("/Users/searleb/Documents/projects/phosphopedia/sp_iso_HUMAN_4.9.2015_UP000005640.fasta");
		File f=new File("/Users/bsearle/Documents/iarpa/2018_1203/10pFDR_DB/IARPA_20181214_10_pepFDR.fasta");
		PrintWriter writer=new PrintWriter("/Users/bsearle/Documents/iarpa/2018_1203/10pFDR_DB/IARPA_20181214_10_pepFDR.csv");
		
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		AminoAcidConstants constants=new AminoAcidConstants();
		
		HashSet<String>[] allPeptides=new HashSet[2];
		for (int i=0; i<allPeptides.length; i++) {
			allPeptides[i]=new HashSet<>();
		}
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		for (FastaEntryInterface entry : entries) {
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 7, 30, 1, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
			for (FastaPeptideEntry pep : peptides) {
				for (int pepCharge : new int[] {2, 3}) {
					String seq=pep.getSequence();
					double pepMass=constants.getMass(seq)+MassConstants.oh2;
					double pepChargedMass=(pepMass+MassConstants.protonMass*pepCharge)/pepCharge;

					if (pepChargedMass>(400)&&pepChargedMass<(1000)) {
						if (seq.indexOf('B')>=0||seq.indexOf('J')>=0||seq.indexOf('O')>=0||seq.indexOf('U')>=0||seq.indexOf('X')>=0||seq.indexOf('Z')>=0||seq.indexOf('*')>=0) {
							continue;
						} else {
							allPeptides[pepCharge-2].add(seq);
						}
					}
				}
			}
		}
		
		writer.println("modified_sequence,collision_energy,precursor_charge");
		for (int i=0; i<allPeptides.length; i++) {
			int charge=i+2;
			for (String string : allPeptides[i]) {
				writer.println(string+","+convertNCE(27f, (byte)charge, (byte)2)+","+(charge));
			}
		}
		writer.close();
	}
	
	// http://proteomicsnews.blogspot.com/2014/06/normalized-collision-energy-calculation.html
	private static float convertNCE(float nce, byte charge, byte defaultCharge) {
		return nce*getChargeFactor(defaultCharge)/getChargeFactor(charge);
	}

	private static float getChargeFactor(byte charge) {
		switch (charge) {
			case 1: return 1.0f;
			case 2: return 0.9f;
			case 3: return 0.85f;
			case 4: return 0.8f;
			case 5: return 0.75f;
			default: return 0.75f;
		}
	}
	
	public static void main4(String[] args) {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		//File f=new File("/Users/searleb/Documents/projects/phosphopedia/sp_iso_HUMAN_4.9.2015_UP000005640.fasta");
		File f=new File("/Users/searleb/Documents/school/uniprot-9606.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		AminoAcidConstants constants=new AminoAcidConstants();
		System.out.println(entries.size());

		int countKR=0;
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		for (FastaEntryInterface entry : entries) {
			int charge=1+getCount(entry.getSequence(), 'K', 'R');
			double mass=constants.getMass(entry.getSequence())+MassConstants.oh2;
			double chargedMass=(mass+MassConstants.protonMass*charge)/charge;
			
			//System.out.println(charge);
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 8, 40, 2, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
			for (FastaPeptideEntry string : peptides) {
				int pepCharge=2;
				double pepMass=constants.getMass(string.getSequence())+MassConstants.oh2;
				double pepChargedMass=(pepMass+MassConstants.protonMass*pepCharge)/pepCharge;
				
				if (pepChargedMass>(665.3204-5)&&pepChargedMass<(665.3204+5)) {
					System.out.println(string);
				}
			}
		}
	}

	public static void main3(String[] args) {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		AminoAcidConstants aaConstants=parameters.getAAConstants();
		//File f=new File("/Users/searleb/Documents/projects/phosphopedia/sp_iso_HUMAN_4.9.2015_UP000005640.fasta");
		File f=new File("/Users/searleb/Documents/school/uniprot-9606.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		Range possiblePrecursors=new Range(400, 1000);
		
		DigestionEnzyme trypsin=DigestionEnzyme.getEnzyme("trypsin");
		DigestionEnzyme gluc=DigestionEnzyme.getEnzyme("glu-c");
		DigestionEnzyme chymotrypsin=DigestionEnzyme.getEnzyme("chymotrypsin");
		DigestionEnzyme lysc=DigestionEnzyme.getEnzyme("lys-c");
		DigestionEnzyme argc=DigestionEnzyme.getEnzyme("arg-c");
		DigestionEnzyme aspn=DigestionEnzyme.getEnzyme("asp-n");
			
		int totalTotalSTY=0;
		int totalVisibleSTY=0;
		int totalTheOnlyOneSTY=0;
		int totalOneOfTwoSTY=0;
		for (FastaEntryInterface entry : entries) {
			int totalSTY=getCount(entry.getSequence(), 'S', 'T', 'Y');
			TIntHashSet visibleSTYIndicies=new TIntHashSet();
			TIntHashSet theOnlyOneSTYIndicies=new TIntHashSet();
			TIntHashSet oneOfTwoSTYIndicies=new TIntHashSet();

			for (DigestionEnzyme enzyme : new DigestionEnzyme[] {trypsin, chymotrypsin}) { //,trypsin,chymotrypsin,gluc,aspn,lysc,argc
				ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 0, 9999999, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
				int peptideIndex=0;
				for (FastaPeptideEntry peptide : peptides) {
					String sequence=peptide.getSequence();
	
					double mass=aaConstants.getMass(sequence)+MassConstants.oh2+79.966331;
					boolean ok=false;
					for (int charge=2; charge<=4; charge++) {
						double chargedMass=(mass+MassConstants.protonMass*charge)/charge;
						if (possiblePrecursors.contains(chargedMass)) {
							ok=true;
							break;
						}
					}
	
					int countSTY=getCount(sequence, 'S', 'T', 'Y');
					
					if (ok) {
						visibleSTYIndicies.addAll(General.add(getAllIndicies(sequence, 'S', 'T', 'Y'), peptideIndex));
						if (countSTY==1) {
							theOnlyOneSTYIndicies.add(peptideIndex+getFirstIndex(sequence, 'S', 'T', 'Y'));
						}
						if (countSTY<=2) {
							oneOfTwoSTYIndicies.addAll(General.add(getAllIndicies(sequence, 'S', 'T', 'Y'), peptideIndex));
						}
					}
					
					peptideIndex+=sequence.length();
				}
			}
			totalTotalSTY+=totalSTY;
			totalVisibleSTY+=visibleSTYIndicies.size();
			totalTheOnlyOneSTY+=theOnlyOneSTYIndicies.size();
			totalOneOfTwoSTY+=oneOfTwoSTYIndicies.size();
		}
		System.out.println(totalTotalSTY+"\t"+totalVisibleSTY+"\t"+totalTheOnlyOneSTY+"\t"+totalOneOfTwoSTY);
	}
	
	public static void main5(String[] args) {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		//File f=new File("/Users/searleb/Documents/projects/phosphopedia/sp_iso_HUMAN_4.9.2015_UP000005640.fasta");
		File f=new File("/Users/searleb/Documents/school/uniprot-9606.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);

		int countBase=0;
		int countNTermProtein=0;
		int countNTermPyroGlu=0;
		int countTryp=0;
		int countMet=0;
		int countSTY=0;
		int countQN=0;
		int countKR=0;
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		for (FastaEntryInterface entry : entries) {
			countNTermProtein++;
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 8, 40, 2, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
			for (FastaPeptideEntry peptide : peptides) {
				String sequence=peptide.getSequence();
				countBase++;
				countNTermProtein++;
				if (sequence.charAt(0)=='Q'||sequence.charAt(0)=='C') {
					countNTermPyroGlu+=2;
				} else {
					countNTermPyroGlu++;
				}
				
				countTryp+=getCombinatorial(sequence, 'W');
				countMet+=getCombinatorial(sequence, 'M');
				countSTY+=getCombinatorial(sequence, 'S', 'T', 'Y');
				countQN+=getCombinatorial(sequence, 'Q', 'N');
				countKR+=getCombinatorial(sequence, 'K', 'R');
			}
		}
		System.out.println(countBase+"\tcountBase");
		System.out.println(countNTermProtein+"\tcountNTermProtein");
		System.out.println(countNTermPyroGlu+"\tcountNTermPyroGlu");
		System.out.println(countTryp+"\tcountTryp");
		System.out.println(countMet+"\tcountMet");
		System.out.println(countQN+"\tcountQN");
		System.out.println(countSTY+"\tcountSTY");
	}
	
	static int getCombinatorial(String sequence, char... target) {
		int num=getCount(sequence, target);
		if (num==0) return 1;
		if (num==1) return 2;
		if (num==2) return 4;
		return (int)(1+num+CombinatoricsUtils.binomialCoefficient(num, 2)+CombinatoricsUtils.binomialCoefficient(num, 3));
	}

	private static int getFirstIndex(String sequence, char... target) {
		int lowestIndex=Integer.MAX_VALUE;
		
		for (int i=0; i<target.length; i++) {
			int index=sequence.indexOf(target[i]);
			if (index>=0&&lowestIndex>index) {
				lowestIndex=index;
			}
		}
		return lowestIndex;
	}

	private static int[] getAllIndicies(String sequence, char... target) {
		TIntArrayList indicies=new TIntArrayList();
		for (int index=0; index<sequence.length(); index++) {
			for (int i=0; i<target.length; i++) {
				if (sequence.charAt(index)==target[i]) {
					indicies.add(index);
				}
			}
		}
		return indicies.toArray();
	}

	private static int getCount(String sequence, char... target) {
		int num=0;
		for (char c : sequence.toCharArray()) {
			for (int i=0; i<target.length; i++) {
				if (c==target[i]) {
					num++;
				}
			}
		}
		return num;
	}
	
	public void testFastaParsing() {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\n"+"ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF\n"
				+"FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV\n"+"ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK\n"
				+"ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR\n"+"RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE\n"
				+"QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV\n"+"LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL\n"
				+"SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV\n"+"AASQAALGL";

		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		assertEquals("ALBU_HUMAN", entry.getAccession());
		assertEquals("MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL", entry.getSequence());

		String ecoli=">gi|16131183|ref|NP_417763.1| 50S ribosomal subunit protein L18 [Escherichia coli str. K-12 substr. MG1655]\n"
				+"MDKKSARIRRATRARRKLQELGATRLVVHRTPRHIYAQVIAPNGSEVLVAASTVEKAIAEQLKYTGNKDA\n"+"AAAVGKAVAERALEKGIKDVSFDRSGFQYHGRVQALADAAREAGLQF\n";

		entry=FastaReader.readFasta(ecoli, "", parameters).get(0);
		assertEquals("MDKKSARIRRATRARRKLQELGATRLVVHRTPRHIYAQVIAPNGSEVLVAASTVEKAIAEQLKYTGNKDAAAAVGKAVAERALEKGIKDVSFDRSGFQYHGRVQALADAAREAGLQF", entry.getSequence());
		ArrayList<FastaPeptideEntry> peptides=DigestionEnzyme.getEnzyme("trypsin").digestProtein(entry, 8, 40, 1, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
		boolean found=false;
		for (FastaPeptideEntry fastaPeptideEntry : peptides) {
			if ("GIKDVSFDR".equals(fastaPeptideEntry.getSequence())) {
				found=true;
				break;
			}
		}
		assertTrue(found);
	}

	public void testFastaReader() throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		InputStream is=getClass().getResourceAsStream("/ecoli-190209-contam_correctNL.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(is, "ecoli-190209-contam_correctNL.fasta", parameters);

		assertEquals(4178, entries.size());
		
		File temp=File.createTempFile("ecoli_", ".fasta"); // test write
		temp.deleteOnExit();
		FastaWriter.writeFasta(temp, entries);
		ArrayList<FastaEntryInterface> writtenEntries=FastaReader.readFasta(temp, parameters);

		assertEquals(4178, writtenEntries.size());
		for (int i=0; i<writtenEntries.size(); i++) {
			assertEquals(entries.get(i).getAccession(), writtenEntries.get(i).getAccession());
			assertEquals(entries.get(i).getSequence(), writtenEntries.get(i).getSequence());
		}
		
		entries=FastaReader.readFasta(temp, "NP_", parameters); // test filter
		assertEquals(3934, entries.size());

		File reverseConcatenated=PercolatorExecutor.getFastaPlusDecoyFile(temp, parameters);
		ArrayList<FastaEntryInterface> reverseConcatenatedEntries=FastaReader.readFasta(reverseConcatenated, parameters);
		assertEquals(4178*2, reverseConcatenatedEntries.size());
		reverseConcatenated.deleteOnExit();

		/*
		TIntObjectHashMap<TFloatArrayList> peptideDefects=new TIntObjectHashMap<TFloatArrayList>();
		
		for (FastaEntry entry : entries) {
			ArrayList<String> peptides=PARAMETERS.getEnzyme().digestProtein(entry.getSequence(), PARAMETERS.getMinPeptideLength(), PARAMETERS.getMaxPeptideLength(), PARAMETERS.getMaxMissedCleavages());
			for (String sequence : peptides) {
				FragmentationModel model=PeptideUtils.getPeptideModel(sequence);
				double[] ions=model.getPrimaryIons(PARAMETERS.getFragType());
				for (double d : ions) {
					d=(d+1.00727646681290)/2;
					int nominalMass=(int)d;
					float defect=(float)(d-nominalMass);
					TFloatArrayList list=peptideDefects.get(nominalMass);
					if (list==null) {
						list=new TFloatArrayList();
						peptideDefects.put(nominalMass, list);
					}
					list.add(defect);
				}
			}
		}
		
		int[] keys=peptideDefects.keys();
		Arrays.sort(keys);
		for (int nominal : keys) {
			float[] defects=peptideDefects.get(nominal).toArray();
			Arrays.sort(defects);
			System.out.println(nominal+"\t"+defects[(int)(defects.length*0.05f)]+"\t"+defects[(int)(defects.length*0.25f)]+"\t"+defects[(int)(defects.length*0.5f)]+"\t"+defects[(int)(defects.length*0.75f)]+"\t"+defects[(int)(defects.length*0.95f)]);
		}
		*/
	}
	
	/*
	 * @MoMo
	 * test reading Peff format and using ExtendedFastaEntry 
	 */
	public void testFastaReaderForPeff() throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		InputStream is=getClass().getResourceAsStream("/nextprot2017_testPEFF1.0rc25_small.peff");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(new BufferedReader(new InputStreamReader(is)), "nextprot2017_testPEFF1.0rc25_small.peff", null, true, parameters);
		assertEquals(25, entries.size());

		for (FastaEntryInterface entry : entries) {
			if (!(entry instanceof ExtendedFastaEntry)) {
				throw new Exception("Error occured when reading peff file, each entry should be an ExtendedFastaEntry object");
			}
		}
	}
}
