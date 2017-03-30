package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;

import org.apache.commons.math3.util.CombinatoricsUtils;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.FastaWriter;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import junit.framework.TestCase;

public class FastaReaderTest extends TestCase {
	
	public static void main(String[] args) {
		File f=new File("/Users/searleb/Downloads/Saccharomyces_cerevisiae_sprot_032417.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f);

		int countBase=0;
		int countNTermProtein=0;
		int countNTermPyroGlu=0;
		int countTryp=0;
		int countMet=0;
		int countSTY=0;
		int countQN=0;
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		for (FastaEntryInterface entry : entries) {
			countNTermProtein++;
			ArrayList<String> peptides=enzyme.digestProtein(entry.getSequence(), 8, 40, 2);
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
		int num=0;
		for (char c : sequence.toCharArray()) {
			for (int i=0; i<target.length; i++) {
				if (c==target[i]) {
					num++;
				}
			}
		}
		if (num==0) return 1;
		if (num==1) return 2;
		if (num==2) return 4;
		return (int)(1+num+CombinatoricsUtils.binomialCoefficient(num, 2)+CombinatoricsUtils.binomialCoefficient(num, 3));
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
		ArrayList<String> peptides=DigestionEnzyme.getEnzyme("trypsin").digestProtein(entry.getSequence(), 8, 40, 1);
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

		/*
		TIntObjectHashMap<TFloatArrayList> peptideDefects=new TIntObjectHashMap<TFloatArrayList>();
		
		for (FastaEntry entry : entries) {
			ArrayList<String> peptides=PARAMETERS.getEnzyme().digestProtein(entry.getSequence(), PARAMETERS.getMinPeptideLength(), PARAMETERS.getMaxPeptideLength(), PARAMETERS.getMaxMissedCleavages());
			for (String sequence : peptides) {
				FragmentationModel model=new FragmentationModel(sequence);
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
}
