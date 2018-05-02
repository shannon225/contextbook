package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.apache.commons.math3.util.CombinatoricsUtils;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.ExtendedFastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filewriters.FastaWriter;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.map.hash.TCharDoubleHashMap;
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
		PecanSearchParameters parameters=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false, true);
		//File f=new File("/Users/searleb/Documents/projects/phosphopedia/sp_iso_HUMAN_4.9.2015_UP000005640.fasta");
		File f=new File("/Users/searleb/Documents/chromatogram_library_manuscript/real_pecan/cerevisiae_orf_trans_all.fasta");
		ArrayList<FastaEntryInterface> targetProteins=FastaReader.readFasta(f);
		
		//PrintWriter writer=new PrintWriter("/Users/searleb/Documents/chromatogram_library_manuscript/sp_iso_HUMAN_4.9.2015_UP000005640.peptides.txt");
		PrintWriter writer=new PrintWriter("/Users/searleb/Documents/chromatogram_library_manuscript/real_pecan/cerevisiae_orf_trans_all.peptides.txt");
		for (FastaEntryInterface entry : targetProteins) {
			ArrayList<String> peptides=parameters.getEnzyme().digestProtein(entry, parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages(), parameters.getAAConstants());
			for (String peptide : peptides) {
				writer.println(entry.getAccession()+"\t"+PeptideUtils.getPeptideSeq(peptide));
			}
		}
		writer.flush();
		writer.close();
	}
	
	public static void main(String[] args) {
		//File f=new File("/Users/searleb/Documents/projects/phosphopedia/sp_iso_HUMAN_4.9.2015_UP000005640.fasta");
		File f=new File("/Users/searleb/Documents/school/uniprot-9606.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f);
		AminoAcidConstants constants=new AminoAcidConstants();
		System.out.println(entries.size());

		int countKR=0;
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		for (FastaEntryInterface entry : entries) {
			int charge=1+getCount(entry.getSequence(), 'K', 'R');
			double mass=constants.getMass(entry.getSequence())+MassConstants.oh2;
			double chargedMass=(mass+MassConstants.protonMass*charge)/charge;
			
			//System.out.println(charge);
			ArrayList<String> peptides=enzyme.digestProtein(entry, 8, 40, 2, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
			for (String string : peptides) {
				int pepCharge=2;
				double pepMass=constants.getMass(string)+MassConstants.oh2;
				double pepChargedMass=(pepMass+MassConstants.protonMass*pepCharge)/pepCharge;
				
				if (pepChargedMass>(665.3204-5)&&pepChargedMass<(665.3204+5)) {
					System.out.println(string);
				}
			}
		}
	}
	
	public static void main3(String[] args) {
		//File f=new File("/Users/searleb/Documents/projects/phosphopedia/sp_iso_HUMAN_4.9.2015_UP000005640.fasta");
		File f=new File("/Users/searleb/Documents/school/uniprot-9606.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f);

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
			ArrayList<String> peptides=enzyme.digestProtein(entry, 8, 40, 2, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
			for (String sequence : peptides) {
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
				countQN+=getCombinatorial(sequence, 'K', 'R');
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
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\n"+"ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF\n"
				+"FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV\n"+"ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK\n"
				+"ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR\n"+"RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE\n"
				+"QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV\n"+"LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL\n"
				+"SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV\n"+"AASQAALGL";

		FastaEntryInterface entry=FastaReader.readFasta(bsa, "").get(0);
		assertEquals("ALBU_HUMAN", entry.getAccession());
		assertEquals("MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL", entry.getSequence());

		String ecoli=">gi|16131183|ref|NP_417763.1| 50S ribosomal subunit protein L18 [Escherichia coli str. K-12 substr. MG1655]\n"
				+"MDKKSARIRRATRARRKLQELGATRLVVHRTPRHIYAQVIAPNGSEVLVAASTVEKAIAEQLKYTGNKDA\n"+"AAAVGKAVAERALEKGIKDVSFDRSGFQYHGRVQALADAAREAGLQF\n";

		entry=FastaReader.readFasta(ecoli, "").get(0);
		assertEquals("MDKKSARIRRATRARRKLQELGATRLVVHRTPRHIYAQVIAPNGSEVLVAASTVEKAIAEQLKYTGNKDAAAAVGKAVAERALEKGIKDVSFDRSGFQYHGRVQALADAAREAGLQF", entry.getSequence());
		ArrayList<String> peptides=DigestionEnzyme.getEnzyme("trypsin").digestProtein(entry, 8, 40, 1, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
		assertTrue(peptides.contains("GIKDVSFDR"));
	}

	public void testFastaReader() throws Exception {
		InputStream is=getClass().getResourceAsStream("/ecoli-190209-contam_correctNL.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(is, "ecoli-190209-contam_correctNL.fasta");

		assertEquals(4178, entries.size());
		
		File temp=File.createTempFile("ecoli_", ".fasta"); // test write
		temp.deleteOnExit();
		FastaWriter.writeFasta(temp, entries);
		ArrayList<FastaEntryInterface> writtenEntries=FastaReader.readFasta(temp);

		assertEquals(4178, writtenEntries.size());
		for (int i=0; i<writtenEntries.size(); i++) {
			assertEquals(entries.get(i).getAccession(), writtenEntries.get(i).getAccession());
			assertEquals(entries.get(i).getSequence(), writtenEntries.get(i).getSequence());
		}
		
		entries=FastaReader.readFasta(temp, "NP_"); // test filter
		assertEquals(3934, entries.size());

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		File reverseConcatenated=PercolatorExecutor.getFastaPlusDecoyFile(temp, parameters);
		ArrayList<FastaEntryInterface> reverseConcatenatedEntries=FastaReader.readFasta(reverseConcatenated);
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
		InputStream is=getClass().getResourceAsStream("/nextprot2017_testPEFF1.0rc25_small.peff");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(new BufferedReader(new InputStreamReader(is)), "nextprot2017_testPEFF1.0rc25_small.peff", null, true);
		assertEquals(25, entries.size());

		for (FastaEntryInterface entry : entries) {
			if (!(entry instanceof ExtendedFastaEntry)) {
				throw new Exception("Error occured when reading peff file, each entry should be an ExtendedFastaEntry object");
			}
		}
	}
}
