package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.InputStream;
import java.util.ArrayList;

import junit.framework.TestCase;

public class FastaReaderTest extends TestCase {
	//private static final SearchParameters PARAMETERS=new SearchParameters(FragmentationType.CID, new MassTolerance(50), new MassTolerance(50), DigestionEnzyme.getEnzyme("trypsin"));
	
	public void testFastaParsing() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\n"+"ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF\n"
				+"FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV\n"+"ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK\n"
				+"ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR\n"+"RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE\n"
				+"QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV\n"+"LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL\n"
				+"SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV\n"+"AASQAALGL";

		FastaEntry entry=FastaReader.readFasta(bsa, "").get(0);
		assertEquals(">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2", entry.getAnnotation());
		assertEquals("MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL", entry.getSequence());
	}

	public void testFastaReader() {
		InputStream is=getClass().getResourceAsStream("/ecoli-190209-contam_correctNL.fasta");
		ArrayList<FastaEntry> entries=FastaReader.readFasta(is, "ecoli-190209-contam_correctNL.fasta");

		assertEquals(4178, entries.size());

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
