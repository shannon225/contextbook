package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.StringUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.procedure.TObjectIntProcedure;
import gnu.trove.set.hash.TCharHashSet;
import junit.framework.TestCase;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.*;

import static edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme.AAs;
import static edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme.getAvailableEnzymes;

public class DigestionEnzymeTest extends TestCase {
	
	public static void main(String[] args) {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		DigestionEnzyme enzyme = DigestionEnzyme.getEnzyme("Trypsin");
		File f=new File("/Users/searleb/Downloads/mus_musculus_reviewed_uniprot.fasta");
		//File f=new File("/Users/searleb/Downloads/homo_sapiens_reviewed_uniprot.fasta");

		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		for (FastaEntryInterface entry : entries) {
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 7, 25, 0, parameters.getAAConstants(), false);
			for (FastaPeptideEntry peptide : peptides) {
				if (peptide.getSequence().indexOf('X')>=0||peptide.getSequence().indexOf('B')>=0||peptide.getSequence().indexOf('Z')>=0||peptide.getSequence().indexOf('U')>=0) {
					continue;
				}
				System.out.println(peptide.getSequence());
			}
		}
	}
	
	public static void mainq(String[] args) {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		File f=new File("/Users/searleb/Downloads/uniprotkb_proteome_UP000005640_2024_01_07.fasta");

		//File f=new File("/Users/searleb/Downloads/uniprotkb_reviewed_true_AND_model_organ_2024_01_13.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);

		TObjectIntHashMap<String> map=new TObjectIntHashMap<String>();
		for (FastaEntryInterface entry : entries) {
			ArrayList<FastaPeptideEntry> peptides=DigestionEnzyme.getEnzyme("Trypsin").digestProtein(entry, 7, 35, 0, parameters.getAAConstants(), false);
			for (FastaPeptideEntry peptide : peptides) {
				String sequence = peptide.getSequence().replace('I', 'L');
				map.put(sequence, map.get(sequence)+1);
			}
		}
		
		final int[] counts=new int[11];
		map.forEachEntry(new TObjectIntProcedure<String>() {
			@Override
			public boolean execute(String a, int b) {
				if (b>=counts.length) b=counts.length-1;
				counts[b]++;
				return true;
			}
		});
		
		for (int i = 0; i < counts.length; i++) {
			System.out.println(i+"\t"+counts[i]);
		}
		System.out.println(1.0f-(counts[1]/(float)General.sum(counts)));
	}
	
	public static void main2(String[] args) {
		String cdk1="MAEAPASPAPLSPLEVELDPEFEPQSRPRSCTWPLQRPELQASPAKPSGETAADSMIPEE" + 
				"EDDEDDEDGGGRAGSAMAIGGGGGSGTLGSGLLLEDSARVLAPGGQDPGSGPATAAGGLS" + 
				"GGTQALLQPQQPLPPPQPGAAGGSGQPRKCSSRRNAWGNLSYADLITRAIESSPDKRLTL" + 
				"SQIYEWMVRCVPYFKDKGDSNSSAGWKNSIRHNLSLHSRFMRVQNEGTGKSSWWIINPDG" + 
				"GKSGKAPRRRAVSMDNSNKYTKSRGRAAKKKAALQTAPESADDSPSQLSKWPGSPTSRSS" + 
				"DELDAWTDFRSRTNSNASTVSGRLSPIMASTELDEVQDDDAPLSPMLYSSSASLSPSVSK" + 
				"PCTVELPRLTDMAGTMNLNDGLTENLMDDLLDNITLPPSQPSPTGGLMQRSSSFPYTTKG" + 
				"SGLGSPTSSFNSTVFGPSSLNSLRQSPMQTIQENKPATFSSMSHYGNQTLQDLLTSDSLS" + 
				"HSDVMMTQSDPLMSQASTAVSAQNSRRNVMLRNDPMMSFAAQPNQGSLVNQNLLHHQHQT" + 
				"QGALGGSRALSNSVSNMGLSESSSLGSAKHQQQSPVSQSMQTLSDSLSGSSLYSTSANLP" + 
				"VMGHEKFPSDLDLDMFNGSLECDMESIIRSELMDADGLDFNFDSLISTQNVVGLNVGNFT" + 
				"GAKQASSQSWVPG";
		FastaEntry entry=new FastaEntry(cdk1);

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		AminoAcidConstants aaConstants=parameters.getAAConstants();
		Range possiblePrecursors=new Range(400, 1000);
		
		DigestionEnzyme trypsin=DigestionEnzyme.getEnzyme("trypsin");
		DigestionEnzyme gluc=DigestionEnzyme.getEnzyme("glu-c");
		DigestionEnzyme chymotrypsin=DigestionEnzyme.getEnzyme("chymotrypsin");
		System.out.println(StringUtils.wrap(cdk1, 100));
		for (DigestionEnzyme enzyme : new DigestionEnzyme[] {trypsin, chymotrypsin, gluc}) { //,trypsin,chymotrypsin,gluc,aspn,lysc,argc
			StringBuilder sb=new StringBuilder();
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 0, 9999999, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
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
				if (ok) {
					sb.append(sequence);
				} else {
					char[] empty=new char[sequence.length()];
					Arrays.fill(empty, ' ');
					sb.append(empty);
				}
			}
			System.out.println(StringUtils.wrap(sb.toString(), 100));
		}
	}
	public void testReverse() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPENECFLQH\n";

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		String sequence=entry.getSequence();
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		AminoAcidConstants aminoAcidConstants = new AminoAcidConstants();
		String reversed=enzyme.reverseProtein(sequence, aminoAcidConstants);
		assertEquals("MSYASSFLFLLSIFTVWRGFVRRDHAKSHAVERFKDFNEEGLKAVHDEFPCQQLYQAFAILVLKLAFETVENVKTDCNEASEDAVCKSDGFLTHLKLLTAVTCREACCDAMEGYTKQQLFCENEPEH", reversed);
	}
	

	public void testNTermMHandling() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPENECFLQH\n";
		
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		HashSet<String> expected=new HashSet<String>();
		expected.add("MWVTFISLLFLFSSAYSR");
		expected.add("WVTFISLLFLFSSAYSR");
		expected.add("DLGEENFK");
		expected.add("ALVLIAFAQYLQQC[57.0214635]PFEDHVK");
		expected.add("LVNEVTEFAK");
		expected.add("TC[57.0214635]VADESAENC[57.0214635]DK");
		expected.add("SLHTLFGDK");
		expected.add("LC[57.0214635]TVATLR");
		expected.add("ETYGEMADC[57.0214635]C[57.0214635]AK");
		expected.add("QEPENEC[57.0214635]FLQH");

		ModificationMassMap variableMods=new ModificationMassMap(); //"C=14.01565");
		AminoAcidConstants constants=new AminoAcidConstants(AminoAcidConstants.getFixedModsMap("C+57 (Carbamidomethyl)"), variableMods);

		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 40, 0, true, constants, false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}
	}
		

	public void testFixedModsNoNTermMHandling() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPENECFLQH\n";

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		String sequence=entry.getSequence();

		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		HashSet<String> expected=new HashSet<String>();
		expected.add("MWVTFISLLFLFSSAYSR");
		expected.add("DLGEENFK");
		expected.add("ALVLIAFAQYLQQC[57.0214635]PFEDHVK");
		expected.add("LVNEVTEFAK");
		expected.add("TC[57.0214635]VADESAENC[57.0214635]DK");
		expected.add("SLHTLFGDK");
		expected.add("LC[57.0214635]TVATLR");
		expected.add("ETYGEMADC[57.0214635]C[57.0214635]AK");
		expected.add("QEPENEC[57.0214635]FLQH");

		ModificationMassMap variableMods=new ModificationMassMap(); //"C=14.01565");
		AminoAcidConstants constants=new AminoAcidConstants(AminoAcidConstants.getFixedModsMap("C+57 (Carbamidomethyl)"), variableMods);

		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 40, 0, false, constants, false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}

		expected=new HashSet<String>();
		expected.add("MWVTFISLLFLFSSAYSR");
		expected.add("DLGEENFK");
		expected.add("ALVLIAFAQYLQQC[57.0214635]PFEDHVK");
		expected.add("LVNEVTEFAK");
		expected.add("TC[57.0214635]VADESAENC[57.0214635]DK");
		expected.add("SLHTLFGDK");
		expected.add("LC[57.0214635]TVATLR");
		expected.add("ETYGEMADC[57.0214635]C[57.0214635]AK");
		expected.add("QEPENEC[57.0214635]FLQH");

		expected.add("ALVLIAFAQYLQQC[14.01565]PFEDHVK");
		expected.add("TC[57.0214635]VADESAENC[14.01565]DK");
		expected.add("TC[14.01565]VADESAENC[57.0214635]DK");
		expected.add("LC[14.01565]TVATLR");
		expected.add("ETYGEMADC[57.0214635]C[14.01565]AK");
		expected.add("ETYGEMADC[14.01565]C[57.0214635]AK");
		expected.add("QEPENEC[14.01565]FLQH");

		variableMods=new ModificationMassMap("C=14.01565");
		constants=new AminoAcidConstants(AminoAcidConstants.getFixedModsMap("C+57 (Carbamidomethyl)"), variableMods);

		sequences=enzyme.digestProtein(entry, 8, 40, 0, false, constants, false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}
	}

	public void testFixedMods() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPENECFLQH\n";

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		String sequence=entry.getSequence();
		
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		HashSet<String> expected=new HashSet<String>();
		expected.add("MWVTFISLLFLFSSAYSR");
		expected.add("WVTFISLLFLFSSAYSR");
		expected.add("DLGEENFK");
		expected.add("ALVLIAFAQYLQQC[57.0214635]PFEDHVK");
		expected.add("LVNEVTEFAK");
		expected.add("TC[57.0214635]VADESAENC[57.0214635]DK");
		expected.add("SLHTLFGDK");
		expected.add("LC[57.0214635]TVATLR");
		expected.add("ETYGEMADC[57.0214635]C[57.0214635]AK");
		expected.add("QEPENEC[57.0214635]FLQH");
		
		ModificationMassMap variableMods=new ModificationMassMap(); //"C=14.01565");
		AminoAcidConstants constants=new AminoAcidConstants(AminoAcidConstants.getFixedModsMap("C+57 (Carbamidomethyl)"), variableMods);
		
		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 40, 0, constants, false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}
		
		expected=new HashSet<String>();
		expected.add("MWVTFISLLFLFSSAYSR");
		expected.add("WVTFISLLFLFSSAYSR");
		expected.add("DLGEENFK");
		expected.add("ALVLIAFAQYLQQC[57.0214635]PFEDHVK");
		expected.add("LVNEVTEFAK");
		expected.add("TC[57.0214635]VADESAENC[57.0214635]DK");
		expected.add("SLHTLFGDK");
		expected.add("LC[57.0214635]TVATLR");
		expected.add("ETYGEMADC[57.0214635]C[57.0214635]AK");
		expected.add("QEPENEC[57.0214635]FLQH");
		
		expected.add("ALVLIAFAQYLQQC[14.01565]PFEDHVK");
		expected.add("TC[57.0214635]VADESAENC[14.01565]DK");
		expected.add("TC[14.01565]VADESAENC[57.0214635]DK");
		expected.add("LC[14.01565]TVATLR");
		expected.add("ETYGEMADC[57.0214635]C[14.01565]AK");
		expected.add("ETYGEMADC[14.01565]C[57.0214635]AK");
		expected.add("QEPENEC[14.01565]FLQH");
		
		variableMods=new ModificationMassMap("C=14.01565");
		constants=new AminoAcidConstants(AminoAcidConstants.getFixedModsMap("C+57 (Carbamidomethyl)"), variableMods);
		
		sequences=enzyme.digestProtein(entry, 8, 40, 0, constants, false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}
	}

	public void testNTermMHandlingRespectsMinLength() {
		String bsa=">TEST_PROTEIN Which I made up\n"+"MAPKANDMAPR\n";

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		String sequence=entry.getSequence();

		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		HashSet<String> expected=new HashSet<String>();
		expected.add("MAPK");
		expected.add("ANDMAPR");


		ArrayList<FastaPeptideEntry> sequences;

		int minLength;

		minLength = 4;
		sequences = enzyme.digestProtein(entry, minLength, 99999999, 0, true,
				new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);

		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}

		minLength = 3;
		expected.add("APK");
		sequences = enzyme.digestProtein(entry, minLength, 99999999, 0, true,
				new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);

		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}
	}

	public void testNTermMModified() {
		String bsa=">TEST_PROTEIN Which I made up\n"+"MAPEPTIDEKANDAPEPTIDER\n";

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		String sequence=entry.getSequence();

		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		HashSet<String> expected=new HashSet<String>();
		expected.add("[42.0]MAPEPTIDEK");
		expected.add("M[16.0]APEPTIDEK");
		expected.add("MAPEPTIDEK");
		expected.add("APEPTIDEK");
		expected.add("ANDAPEPTIDER");


		ModificationMassMap variableMods=new ModificationMassMap("M=16,a=42");
		AminoAcidConstants constants=new AminoAcidConstants(new TCharDoubleHashMap(), variableMods);

		ArrayList<FastaPeptideEntry> sequences;

		sequences = enzyme.digestProtein(entry, 1, 99999999, 0, true,
				constants, false);

		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}
	}

	public void testDigestProteinNoNTermM() {
		String bsa=">TEST_PROTEIN Which I made up\n"+"APEPTIDEKANDAPEPTIDER\n";

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		String sequence=entry.getSequence();

		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		HashSet<String> expected=new HashSet<String>();
		expected.add("APEPTIDEK");
		expected.add("ANDAPEPTIDER");


		ArrayList<FastaPeptideEntry> sequences;


		sequences = enzyme.digestProtein(entry, 8, 99999999, 0, false,
				new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);

		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}


		sequences = enzyme.digestProtein(entry, 8, 99999999, 0, true,
				new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);

		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}
	}

	public void testModifications() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPENECFLQH\n";

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		String sequence=entry.getSequence();
		
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		HashSet<String> expected=new HashSet<String>();
		expected.add("MWVTFISLLFLFSSAYSR");
		expected.add("WVTFISLLFLFSSAYSR"); // handling n-term methionine
		expected.add("DLGEENFK");
		expected.add("ALVLIAFAQYLQQCPFEDHVK");
		expected.add("LVNEVTEFAK");
		expected.add("TCVADESAENCDK");
		expected.add("SLHTLFGDK");
		expected.add("LCTVATLR");
		expected.add("ETYGEMADCCAK");
		expected.add("MWVTFISLLFLFSSAYSR");
		expected.add("MW[15.994915]VTFISLLFLFSSAYSR");
		expected.add("W[15.994915]VTFISLLFLFSSAYSR"); // handling n-term methionine
		expected.add("[42.010565]MWVTFISLLFLFSSAYSR");
		expected.add("DLGEENFK");
		expected.add("DLGEENFK[8.014199]");
		expected.add("ALVLIAFAQYLQQCPFEDHVK");
		expected.add("ALVLIAFAQYLQQCPFEDHVK[8.014199]");
		expected.add("LVNEVTEFAK");
		expected.add("LVNEVTEFAK[8.014199]");
		expected.add("TCVADESAENCDK");
		expected.add("TCVADESAENCDK[8.014199]");
		expected.add("SLHTLFGDK");
		expected.add("SLHTLFGDK[8.014199]");
		expected.add("LCTVATLR");
		expected.add("ETYGEMADCCAK");
		expected.add("[-18.026549]ETYGEMADCCAK");
		expected.add("ETYGEMADCCAK[8.014199]");
		expected.add("QEPENECFLQH");
		expected.add("[-17.026549]QEPENECFLQH");
		expected.add("QEPENECFLQH[14.0]");
		
		ModificationMassMap variableMods=new ModificationMassMap("K=8.014199,W=15.994915,nE=-18.026549,nQ=-17.026549,nC=-17.026549,a=42.010565,z=14.0");
		AminoAcidConstants constants=new AminoAcidConstants(new TCharDoubleHashMap(), variableMods);
		
		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 40, 0, constants, false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}
	}



	public void testModificationsNoNTermMHandling() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPENECFLQH\n";

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		String sequence=entry.getSequence();

		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		HashSet<String> expected=new HashSet<String>();
		expected.add("MWVTFISLLFLFSSAYSR");
		expected.add("DLGEENFK");
		expected.add("ALVLIAFAQYLQQCPFEDHVK");
		expected.add("LVNEVTEFAK");
		expected.add("TCVADESAENCDK");
		expected.add("SLHTLFGDK");
		expected.add("LCTVATLR");
		expected.add("ETYGEMADCCAK");
		expected.add("MWVTFISLLFLFSSAYSR");
		expected.add("MW[15.994915]VTFISLLFLFSSAYSR");
		expected.add("[42.010565]MWVTFISLLFLFSSAYSR");
		expected.add("DLGEENFK");
		expected.add("DLGEENFK[8.014199]");
		expected.add("ALVLIAFAQYLQQCPFEDHVK");
		expected.add("ALVLIAFAQYLQQCPFEDHVK[8.014199]");
		expected.add("LVNEVTEFAK");
		expected.add("LVNEVTEFAK[8.014199]");
		expected.add("TCVADESAENCDK");
		expected.add("TCVADESAENCDK[8.014199]");
		expected.add("SLHTLFGDK");
		expected.add("SLHTLFGDK[8.014199]");
		expected.add("LCTVATLR");
		expected.add("ETYGEMADCCAK");
		expected.add("[-18.026549]ETYGEMADCCAK");
		expected.add("ETYGEMADCCAK[8.014199]");
		expected.add("QEPENECFLQH");
		expected.add("[-17.026549]QEPENECFLQH");
		expected.add("QEPENECFLQH[14.0]");

		ModificationMassMap variableMods=new ModificationMassMap("K=8.014199,W=15.994915,nE=-18.026549,nQ=-17.026549,nC=-17.026549,a=42.010565,z=14.0");
		AminoAcidConstants constants=new AminoAcidConstants(new TCharDoubleHashMap(), variableMods);

		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 40, 0, false, constants, false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}
	}
	
	public void testXanderCleavages() {
		System.out.println("BOLA2T:");
		String sequence="MELSAEYLREKLQRDLEAEHVEVEDTTLNRCSCSFRVLVVSAKFEGKPLLQRHRFCTE";
		for (DigestionEnzyme enzyme : DigestionEnzyme.getAvailableEnzymes()) {
			ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(new FastaEntry(sequence), 6, 40, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
			for (FastaPeptideEntry string : sequences) {
				if (string.getSequence().indexOf("RHRF")>=0) {
					System.out.println(enzyme.getName()+": "+string);
				}
			}
		}
		
		System.out.println("\nBOLA2F:");
		sequence="MELSAEYLREKLQRDLEAEHVEVEDTTLNRCSCSFRVLVVSAKFEGKPLLQRHSLDPSMTIHCDMVITYGLDQLENCQTCGTDYIISVLNLLTLIVEQINTKLPSSFVEKLFIPSSKLLFLRYHKDKEVVAVAHAVYQAMLSLKNIPVLETAYKLILGEMTCALNNLLHSLQLPEACSEIKHEAFKNHVFNVDNAKFVVKFDLSALTTIGNAKNSSL";
		for (DigestionEnzyme enzyme : DigestionEnzyme.getAvailableEnzymes()) {
			ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(new FastaEntry(sequence), 6, 40, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
			for (FastaPeptideEntry string : sequences) {
				if (string.getSequence().indexOf("RHSL")>=0) {
					System.out.println(enzyme.getName()+": "+string);
				}
			}
		}
	}

	public void testNoEnzymeNoNTermMHandling() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\n"+"ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF\n"
				+"FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV\n"+"ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK\n"
				+"ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR\n"+"RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE\n"
				+"QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV\n"+"LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL\n"
				+"SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV\n"+"AASQAALGL";

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("no enzyme");

		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 99999999, 1, false,
				new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
		assertEquals(1, sequences.size());
		assertEquals(entry.getSequence(), sequences.get(0).getSequence());
	}


	public void testNoEnzyme() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\n"+"ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF\n"
				+"FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV\n"+"ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK\n"
				+"ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR\n"+"RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE\n"
				+"QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV\n"+"LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL\n"
				+"SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV\n"+"AASQAALGL";

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);

		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("no enzyme");

		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 99999999, 1,
				new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
		assertEquals(2, sequences.size());
		assertEquals(entry.getSequence(), sequences.get(0).getSequence());
	}
	
	public void testCountMissedCleavages() {
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		assertEquals(1, enzyme.getNumMissedCleavages("DLGEENFKALVLIAFAQYLQQCPFEDHVK"));
		assertEquals(0, enzyme.getNumMissedCleavages("LVNEVTEFAK"));
		assertEquals(1, enzyme.getNumMissedCleavages("ALVLIAFAQYLQQCPFEDHVKLVNEVTEFAK"));
		assertEquals(2, enzyme.getNumMissedCleavages("LVNEVTEFAKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAK"));
	}
	
	public void testMissedCleavages() {

		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\n"+"ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF\n"
				+"FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV\n"+"ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK\n"
				+"ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR\n"+"RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE\n"
				+"QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV\n"+"LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL\n"
				+"SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV\n"+"AASQAALGL";

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		String sequence=entry.getSequence();
		
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		assertFalse(enzyme.isCutSite('P', 'R'));
		assertFalse(enzyme.isCutSite('Q', 'R'));
		
		assertFalse(enzyme.isCutSite('R', 'P'));
		assertTrue(enzyme.isCutSite('R', 'Q'));

		assertFalse(enzyme.isTargetPreSite('P'));
		assertFalse(enzyme.isTargetPostSite('P'));
		
		assertTrue(enzyme.isTargetPreSite('R'));
		assertFalse(enzyme.isTargetPostSite('R'));
		
		
		HashSet<String> expected=new HashSet<String>();
		expected.add("WVTFISLLFLFSSAYSR");
		expected.add("MKWVTFISLLFLFSSAYSR");
		expected.add("KWVTFISLLFLFSSAYSR"); // handling n-term methionine
		expected.add("WVTFISLLFLFSSAYSRGVFR");
		expected.add("DAHKSEVAHR");
		expected.add("SEVAHRFK");
		expected.add("DLGEENFK");
		expected.add("FKDLGEENFK");
		expected.add("ALVLIAFAQYLQQCPFEDHVK");
		expected.add("DLGEENFKALVLIAFAQYLQQCPFEDHVK");
		expected.add("LVNEVTEFAK");
		expected.add("ALVLIAFAQYLQQCPFEDHVKLVNEVTEFAK");
		expected.add("TCVADESAENCDK");
		expected.add("LVNEVTEFAKTCVADESAENCDK");
		expected.add("SLHTLFGDK");
		expected.add("TCVADESAENCDKSLHTLFGDK");
		expected.add("LCTVATLR");
		expected.add("SLHTLFGDKLCTVATLR");
		expected.add("ETYGEMADCCAK");
		expected.add("LCTVATLRETYGEMADCCAK");
		expected.add("ETYGEMADCCAKQEPER");
		expected.add("NECFLQHK");
		expected.add("QEPERNECFLQHK");
		expected.add("DDNPNLPR");
		expected.add("NECFLQHKDDNPNLPR");
		expected.add("LVRPEVDVMCTAFHDNEETFLK");
		expected.add("DDNPNLPRLVRPEVDVMCTAFHDNEETFLK");
		expected.add("LVRPEVDVMCTAFHDNEETFLKK");
		expected.add("KYLYEIAR");
		expected.add("YLYEIARR");
		expected.add("HPYFYAPELLFFAK");
		expected.add("RHPYFYAPELLFFAK");
		expected.add("HPYFYAPELLFFAKR");
		expected.add("AAFTECCQAADK");
		expected.add("YKAAFTECCQAADK");
		expected.add("AAFTECCQAADKAACLLPK");
		expected.add("AACLLPKLDELR");
		expected.add("LDELRDEGK");
		expected.add("DEGKASSAK");
		expected.add("LKCASLQK");
		expected.add("CASLQKFGER");
		expected.add("AFKAWAVAR");
		expected.add("AWAVARLSQR");
		expected.add("AEFAEVSK");
		expected.add("FPKAEFAEVSK");
		expected.add("AEFAEVSKLVTDLTK");
		expected.add("VHTECCHGDLLECADDR");
		expected.add("LVTDLTKVHTECCHGDLLECADDR");
		expected.add("VHTECCHGDLLECADDRADLAK");
		expected.add("YICENQDSISSK");
		expected.add("ADLAKYICENQDSISSK");
		expected.add("YICENQDSISSKLK");
		expected.add("ECCEKPLLEK");
		expected.add("LKECCEKPLLEK");
		expected.add("SHCIAEVENDEMPADLPSLAADFVESK");
		expected.add("ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESK");
		expected.add("SHCIAEVENDEMPADLPSLAADFVESKDVCK");
		expected.add("DVCKNYAEAK");
		expected.add("DVFLGMFLYEYAR");
		expected.add("NYAEAKDVFLGMFLYEYAR");
		expected.add("DVFLGMFLYEYARR");
		expected.add("HPDYSVVLLLR");
		expected.add("RHPDYSVVLLLR");
		expected.add("HPDYSVVLLLRLAK");
		expected.add("TYETTLEK");
		expected.add("LAKTYETTLEK");
		expected.add("CCAAADPHECYAK");
		expected.add("TYETTLEKCCAAADPHECYAK");
		expected.add("VFDEFKPLVEEPQNLIK");
		expected.add("CCAAADPHECYAKVFDEFKPLVEEPQNLIK");
		expected.add("QNCELFEQLGEYK");
		expected.add("VFDEFKPLVEEPQNLIKQNCELFEQLGEYK");
		expected.add("FQNALLVR");
		expected.add("QNCELFEQLGEYKFQNALLVR");
		expected.add("FQNALLVRYTK");
		expected.add("VPQVSTPTLVEVSR");
		expected.add("KVPQVSTPTLVEVSR");
		expected.add("VPQVSTPTLVEVSRNLGK");
		expected.add("NLGKVGSK");
		expected.add("CCKHPEAK");
		expected.add("MPCAEDYLSVVLNQLCVLHEK");
		expected.add("RMPCAEDYLSVVLNQLCVLHEK");
		expected.add("MPCAEDYLSVVLNQLCVLHEKTPVSDR");
		expected.add("TPVSDRVTK");
		expected.add("CCTESLVNR");
		expected.add("VTKCCTESLVNR");
		expected.add("RPCFSALEVDETYVPK");
		expected.add("CCTESLVNRRPCFSALEVDETYVPK");
		expected.add("EFNAETFTFHADICTLSEK");
		expected.add("RPCFSALEVDETYVPKEFNAETFTFHADICTLSEK");
		expected.add("EFNAETFTFHADICTLSEKER");
		expected.add("QTALVELVK");
		expected.add("KQTALVELVK");
		expected.add("QTALVELVKHKPK");
		expected.add("AVMDDFAAFVEK");
		expected.add("EQLKAVMDDFAAFVEK");
		expected.add("AVMDDFAAFVEKCCK");
		expected.add("ETCFAEEGK");
		expected.add("ADDKETCFAEEGK");
		expected.add("ETCFAEEGKK");
		expected.add("LVAASQAALGL");
		expected.add("KLVAASQAALGL");

		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 40, 1, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}
		
	}

	public void testMissedCleavagesNoNTermMHandling() {

		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\n"+"ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF\n"
				+"FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV\n"+"ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK\n"
				+"ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR\n"+"RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE\n"
				+"QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV\n"+"LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL\n"
				+"SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV\n"+"AASQAALGL";

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		String sequence=entry.getSequence();

		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		assertFalse(enzyme.isCutSite('P', 'R'));
		assertFalse(enzyme.isCutSite('Q', 'R'));

		assertFalse(enzyme.isCutSite('R', 'P'));
		assertTrue(enzyme.isCutSite('R', 'Q'));

		assertFalse(enzyme.isTargetPreSite('P'));
		assertFalse(enzyme.isTargetPostSite('P'));

		assertTrue(enzyme.isTargetPreSite('R'));
		assertFalse(enzyme.isTargetPostSite('R'));


		HashSet<String> expected=new HashSet<String>();
		expected.add("WVTFISLLFLFSSAYSR");
		expected.add("MKWVTFISLLFLFSSAYSR");
		expected.add("WVTFISLLFLFSSAYSRGVFR");
		expected.add("DAHKSEVAHR");
		expected.add("SEVAHRFK");
		expected.add("DLGEENFK");
		expected.add("FKDLGEENFK");
		expected.add("ALVLIAFAQYLQQCPFEDHVK");
		expected.add("DLGEENFKALVLIAFAQYLQQCPFEDHVK");
		expected.add("LVNEVTEFAK");
		expected.add("ALVLIAFAQYLQQCPFEDHVKLVNEVTEFAK");
		expected.add("TCVADESAENCDK");
		expected.add("LVNEVTEFAKTCVADESAENCDK");
		expected.add("SLHTLFGDK");
		expected.add("TCVADESAENCDKSLHTLFGDK");
		expected.add("LCTVATLR");
		expected.add("SLHTLFGDKLCTVATLR");
		expected.add("ETYGEMADCCAK");
		expected.add("LCTVATLRETYGEMADCCAK");
		expected.add("ETYGEMADCCAKQEPER");
		expected.add("NECFLQHK");
		expected.add("QEPERNECFLQHK");
		expected.add("DDNPNLPR");
		expected.add("NECFLQHKDDNPNLPR");
		expected.add("LVRPEVDVMCTAFHDNEETFLK");
		expected.add("DDNPNLPRLVRPEVDVMCTAFHDNEETFLK");
		expected.add("LVRPEVDVMCTAFHDNEETFLKK");
		expected.add("KYLYEIAR");
		expected.add("YLYEIARR");
		expected.add("HPYFYAPELLFFAK");
		expected.add("RHPYFYAPELLFFAK");
		expected.add("HPYFYAPELLFFAKR");
		expected.add("AAFTECCQAADK");
		expected.add("YKAAFTECCQAADK");
		expected.add("AAFTECCQAADKAACLLPK");
		expected.add("AACLLPKLDELR");
		expected.add("LDELRDEGK");
		expected.add("DEGKASSAK");
		expected.add("LKCASLQK");
		expected.add("CASLQKFGER");
		expected.add("AFKAWAVAR");
		expected.add("AWAVARLSQR");
		expected.add("AEFAEVSK");
		expected.add("FPKAEFAEVSK");
		expected.add("AEFAEVSKLVTDLTK");
		expected.add("VHTECCHGDLLECADDR");
		expected.add("LVTDLTKVHTECCHGDLLECADDR");
		expected.add("VHTECCHGDLLECADDRADLAK");
		expected.add("YICENQDSISSK");
		expected.add("ADLAKYICENQDSISSK");
		expected.add("YICENQDSISSKLK");
		expected.add("ECCEKPLLEK");
		expected.add("LKECCEKPLLEK");
		expected.add("SHCIAEVENDEMPADLPSLAADFVESK");
		expected.add("ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESK");
		expected.add("SHCIAEVENDEMPADLPSLAADFVESKDVCK");
		expected.add("DVCKNYAEAK");
		expected.add("DVFLGMFLYEYAR");
		expected.add("NYAEAKDVFLGMFLYEYAR");
		expected.add("DVFLGMFLYEYARR");
		expected.add("HPDYSVVLLLR");
		expected.add("RHPDYSVVLLLR");
		expected.add("HPDYSVVLLLRLAK");
		expected.add("TYETTLEK");
		expected.add("LAKTYETTLEK");
		expected.add("CCAAADPHECYAK");
		expected.add("TYETTLEKCCAAADPHECYAK");
		expected.add("VFDEFKPLVEEPQNLIK");
		expected.add("CCAAADPHECYAKVFDEFKPLVEEPQNLIK");
		expected.add("QNCELFEQLGEYK");
		expected.add("VFDEFKPLVEEPQNLIKQNCELFEQLGEYK");
		expected.add("FQNALLVR");
		expected.add("QNCELFEQLGEYKFQNALLVR");
		expected.add("FQNALLVRYTK");
		expected.add("VPQVSTPTLVEVSR");
		expected.add("KVPQVSTPTLVEVSR");
		expected.add("VPQVSTPTLVEVSRNLGK");
		expected.add("NLGKVGSK");
		expected.add("CCKHPEAK");
		expected.add("MPCAEDYLSVVLNQLCVLHEK");
		expected.add("RMPCAEDYLSVVLNQLCVLHEK");
		expected.add("MPCAEDYLSVVLNQLCVLHEKTPVSDR");
		expected.add("TPVSDRVTK");
		expected.add("CCTESLVNR");
		expected.add("VTKCCTESLVNR");
		expected.add("RPCFSALEVDETYVPK");
		expected.add("CCTESLVNRRPCFSALEVDETYVPK");
		expected.add("EFNAETFTFHADICTLSEK");
		expected.add("RPCFSALEVDETYVPKEFNAETFTFHADICTLSEK");
		expected.add("EFNAETFTFHADICTLSEKER");
		expected.add("QTALVELVK");
		expected.add("KQTALVELVK");
		expected.add("QTALVELVKHKPK");
		expected.add("AVMDDFAAFVEK");
		expected.add("EQLKAVMDDFAAFVEK");
		expected.add("AVMDDFAAFVEKCCK");
		expected.add("ETCFAEEGK");
		expected.add("ADDKETCFAEEGK");
		expected.add("ETCFAEEGKK");
		expected.add("LVAASQAALGL");
		expected.add("KLVAASQAALGL");

		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 40, 1, false,
				new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}

	}


	public void testEnzymes() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\n"+"ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF\n"
				+"FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV\n"+"ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK\n"
				+"ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR\n"+"RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE\n"
				+"QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV\n"+"LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL\n"
				+"SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV\n"+"AASQAALGL";

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		String sequence=entry.getSequence();
		
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		assertFalse(enzyme.isCutSite('P', 'R'));
		assertFalse(enzyme.isCutSite('Q', 'R'));
		
		assertFalse(enzyme.isCutSite('R', 'P'));
		assertTrue(enzyme.isCutSite('R', 'Q'));

		assertFalse(enzyme.isTargetPreSite('P'));
		assertFalse(enzyme.isTargetPostSite('P'));
		
		assertTrue(enzyme.isTargetPreSite('R'));
		assertFalse(enzyme.isTargetPostSite('R'));
		
		
		HashSet<String> expected=new HashSet<String>();
		expected.add("WVTFISLLFLFSSAYSR");
		expected.add("DLGEENFK");
		expected.add("ALVLIAFAQYLQQCPFEDHVK");
		expected.add("LVNEVTEFAK");
		expected.add("TCVADESAENCDK");
		expected.add("SLHTLFGDK");
		expected.add("LCTVATLR");
		expected.add("ETYGEMADCCAK");
		expected.add("NECFLQHK");
		expected.add("DDNPNLPR");
		expected.add("LVRPEVDVMCTAFHDNEETFLK");
		expected.add("HPYFYAPELLFFAK");
		expected.add("AAFTECCQAADK");
		expected.add("AEFAEVSK");
		expected.add("VHTECCHGDLLECADDR");
		expected.add("YICENQDSISSK");
		expected.add("ECCEKPLLEK");
		expected.add("SHCIAEVENDEMPADLPSLAADFVESK");
		expected.add("DVFLGMFLYEYAR");
		expected.add("HPDYSVVLLLR");
		expected.add("TYETTLEK");
		expected.add("CCAAADPHECYAK");
		expected.add("VFDEFKPLVEEPQNLIK");
		expected.add("QNCELFEQLGEYK");
		expected.add("FQNALLVR");
		expected.add("VPQVSTPTLVEVSR");
		expected.add("MPCAEDYLSVVLNQLCVLHEK");
		expected.add("CCTESLVNR");
		expected.add("RPCFSALEVDETYVPK");
		expected.add("EFNAETFTFHADICTLSEK");
		expected.add("QTALVELVK");
		expected.add("AVMDDFAAFVEK");
		expected.add("ETCFAEEGK");
		expected.add("LVAASQAALGL");
		
		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 40, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}
		expected.clear();
		
		enzyme=DigestionEnzyme.getEnzyme("Chymotrypsin");
		expected.add("RRDAHKSEVAHRF");
		expected.add("KDLGEENF");
		expected.add("KALVLIAF");
		expected.add("EDHVKLVNEVTEF");
		expected.add("AKTCVADESAENCDKSLHTLF");
		expected.add("GDKLCTVATLRETY");
		expected.add("GEMADCCAKQEPERNECF");
		expected.add("LQHKDDNPNLPRLVRPEVDVMCTAF");
		expected.add("EIARRHPY");
		expected.add("AVARLSQRFPKAEF");
		expected.add("AEVSKLVTDLTKVHTECCHGDLLECADDRADLAKY");
		expected.add("VESKDVCKNY");
		expected.add("SVVLLLRLAKTY");
		expected.add("ETTLEKCCAAADPHECY");
		expected.add("KPLVEEPQNLIKQNCELF");
		expected.add("QNALLVRY");
		expected.add("LSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCF");
		expected.add("SALEVDETY");
		expected.add("VEKCCKADDKETCF");
		expected.add("AEEGKKLVAASQAALGL");

		sequences=enzyme.digestProtein(entry, 8, 40, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
		//assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			System.out.println(peptide.getSequence());
			//assertTrue(expected.contains(peptide.getSequence()));
		}
		expected.clear();
		
		enzyme=DigestionEnzyme.getEnzyme("lys-n");
		expected.add("KWVTFISLLFLFSSAYSRGVFRRDAH");
		expected.add("KSEVAHRF");
		expected.add("KDLGEENF");
		expected.add("KALVLIAFAQYLQQCPFEDHV");
		expected.add("KLVNEVTEFA");
		expected.add("KTCVADESAENCD");
		expected.add("KSLHTLFGD");
		expected.add("KLCTVATLRETYGEMADCCA");
		expected.add("KQEPERNECFLQH");
		expected.add("KDDNPNLPRLVRPEVDVMCTAFHDNEETFL");
		expected.add("KYLYEIARRHPYFYAPELLFFA");
		expected.add("KAAFTECCQAAD");
		expected.add("KLDELRDEG");
		expected.add("KAWAVARLSQRFP");
		expected.add("KAEFAEVS");
		expected.add("KVHTECCHGDLLECADDRADLA");
		expected.add("KYICENQDSISS");
		expected.add("KSHCIAEVENDEMPADLPSLAADFVES");
		expected.add("KDVFLGMFLYEYARRHPDYSVVLLLRLA");
		expected.add("KTYETTLE");
		expected.add("KCCAAADPHECYA");
		expected.add("KPLVEEPQNLI");
		expected.add("KQNCELFEQLGEY");
		expected.add("KFQNALLVRYT");
		expected.add("KVPQVSTPTLVEVSRNLG");
		expected.add("KRMPCAEDYLSVVLNQLCVLHE");
		expected.add("KTPVSDRVT");
		expected.add("KCCTESLVNRRPCFSALEVDETYVP");
		expected.add("KEFNAETFTFHADICTLSE");
		expected.add("KQTALVELV");
		expected.add("KAVMDDFAAFVE");
		expected.add("KETCFAEEG");
		expected.add("KLVAASQAALGL");


		sequences=enzyme.digestProtein(entry, 8, 40, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry peptide : sequences) {
			assertTrue(expected.contains(peptide.getSequence()));
		}
	}

	//@MoMo 
	public void testDigestionWithVariantsNoNTermMHandling() {
		String fakeTTR=">nxp:NX_P02766-1 \\DbUniqueId=NX_P02766-1 \\PName=Transthyretin isoform Iso 1 \\GName=TTR \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=147 \\SV=266 \\EV=656 \\PE=1 \\ModResPsi=(62|MOD:00041|L-gamma-carboxyglutamic acid)(69|MOD:00047|O-phospho-L-threonine)(72|MOD:00046|O-phospho-L-serine) \\ModRes=(118||N-linked (GlcNAc...)) "
				+"\\VariantSimple=(2|G)(5|H)(5|C)(8|H)(112|*)(137|T)(146|R)(147|D)(147|*)(56|P) "
				+"\\VariantComplex=(70|75|)(110|110|DK)(110|110|APT)(142|142|)(20|40|) \\Processed=(1|20|signal peptide)(21|147|mature protein)"
				+"\nMASHRLLLLCLAGLVFVSEAGPTGTGESKCPLMVKVLDAVRGSPAINVAVHVFRKAADDTWEPFASGKTSESGELHGLTTEEEFVEGIYKVEIDTKSYWKALGISPFHEHAEVVFTANDSGPRRYTIAALLSPYSYSTTAVVTNPKE";
		HashSet<String> expected=new HashSet<String>();
		// peptides from standard sequence
		expected.add("LLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("GSPAINVAVHVFR");
		expected.add("AADDTWEPFASGK");
		expected.add("TSESGELHGLTTEEEFVEGIYK");
		expected.add("ALGISPFHEHAEVVFTANDSGPR");
		expected.add("YTIAALLSPYSYSTTAVVTNPK");

		// simple variant
		expected.add("MASHCLLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("MASHHLLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("LLHLCLAGLVFVSEAGPTGTGESK");
		expected.add("YTIAALLSPYSYTTTAVVTNPK");
		expected.add("YTIAALLSPYSYSTTAVVTNPR");
		expected.add("ALGISPFHEHA");
		expected.add("KPADDTWEPFASGK");

		// complex variant
		expected.add("THGLTTEEEFVEGIYK");
		expected.add("ALGISPFHEDK");
		expected.add("AEVVFTANDSGPR");
		expected.add("ALGISPFHEAPTAEVVFTANDSGPR");
		expected.add("YTIAALLSPYSYSTTAVTNPK");
		expected.add("LLLLCLAGLVFVSER");

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(fakeTTR.getBytes(StandardCharsets.UTF_8)))), "", "", true, parameters).get(0);

		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");

		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 40, 0, false,
				new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry pep : sequences) {
			if (!expected.contains(pep.getSequence())) {
				System.out.println("generated but not expected: "+pep.getSequence());
			}
			expected.remove(pep.getSequence());
		}

		for (String pep : expected) {
			System.out.println("expected but not generated: "+pep);
		}
		assertEquals(0, expected.size());
	}

	//@MoMo
	public void testDigestionWithVariants() {
		String fakeTTR=">nxp:NX_P02766-1 \\DbUniqueId=NX_P02766-1 \\PName=Transthyretin isoform Iso 1 \\GName=TTR \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=147 \\SV=266 \\EV=656 \\PE=1 \\ModResPsi=(62|MOD:00041|L-gamma-carboxyglutamic acid)(69|MOD:00047|O-phospho-L-threonine)(72|MOD:00046|O-phospho-L-serine) \\ModRes=(118||N-linked (GlcNAc...)) "
				+"\\VariantSimple=(2|G)(5|H)(5|C)(8|H)(112|*)(137|T)(146|R)(147|D)(147|*)(56|P) "
				+"\\VariantComplex=(70|75|)(110|110|DK)(110|110|APT)(142|142|)(20|40|) \\Processed=(1|20|signal peptide)(21|147|mature protein)"
				+"\nMASHRLLLLCLAGLVFVSEAGPTGTGESKCPLMVKVLDAVRGSPAINVAVHVFRKAADDTWEPFASGKTSESGELHGLTTEEEFVEGIYKVEIDTKSYWKALGISPFHEHAEVVFTANDSGPRRYTIAALLSPYSYSTTAVVTNPKE";
		HashSet<String> expected=new HashSet<String>();
		// peptides from standard sequence
		expected.add("LLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("GSPAINVAVHVFR");
		expected.add("AADDTWEPFASGK");
		expected.add("TSESGELHGLTTEEEFVEGIYK");
		expected.add("ALGISPFHEHAEVVFTANDSGPR");
		expected.add("YTIAALLSPYSYSTTAVVTNPK");

		// simple variant
		expected.add("MASHCLLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("MASHHLLLLCLAGLVFVSEAGPTGTGESK");

		expected.add("ASHCLLLLCLAGLVFVSEAGPTGTGESK");// handle N-term M
		expected.add("ASHHLLLLCLAGLVFVSEAGPTGTGESK");// handle N-term M

		expected.add("LLHLCLAGLVFVSEAGPTGTGESK");
		expected.add("YTIAALLSPYSYTTTAVVTNPK");
		expected.add("YTIAALLSPYSYSTTAVVTNPR");
		expected.add("ALGISPFHEHA");
		expected.add("KPADDTWEPFASGK");

		// complex variant
		expected.add("THGLTTEEEFVEGIYK");
		expected.add("ALGISPFHEDK");
		expected.add("AEVVFTANDSGPR");
		expected.add("ALGISPFHEAPTAEVVFTANDSGPR");
		expected.add("YTIAALLSPYSYSTTAVTNPK");
		expected.add("LLLLCLAGLVFVSER");

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(fakeTTR.getBytes(StandardCharsets.UTF_8)))), "", "", true, parameters).get(0);

		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");

		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 40, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry pep : sequences) {
			if (!expected.contains(pep.getSequence())) {
				System.out.println("generated but not expected: "+pep.getSequence());
			}
			expected.remove(pep.getSequence());
		}

		for (String pep : expected) {
			System.out.println("expected but not generated: "+pep);
		}
		assertEquals(0, expected.size());
	}

	//@MoMo 
	public void testMissedCleavageWithVariants() {
		String fakeTTR=">nxp:NX_P02766-1 \\DbUniqueId=NX_P02766-1 \\PName=Transthyretin isoform Iso 1 \\GName=TTR \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=147 \\SV=266 \\EV=656 \\PE=1 \\ModResPsi=(62|MOD:00041|L-gamma-carboxyglutamic acid)(69|MOD:00047|O-phospho-L-threonine)(72|MOD:00046|O-phospho-L-serine) \\ModRes=(118||N-linked (GlcNAc...)) "
				+"\\VariantSimple=(2|G)(5|H)(5|C)(8|H)(112|*)(137|T)(146|R)(147|D)(147|*)(56|P) "
				+"\\VariantComplex=(70|75|)(110|110|DK)(110|110|APT)(142|142|)(20|40|)(135|147|) \\Processed=(1|20|signal peptide)(21|147|mature protein)"
				+"\nMASHRLLLLCLAGLVFVSEAGPTGTGESKCPLMVKVLDAVRGSPAINVAVHVFRKAADDTWEPFASGKTSESGELHGLTTEEEFVEGIYKVEIDTKSYWKALGISPFHEHAEVVFTANDSGPRRYTIAALLSPYSYSTTAVVTNPKE";
		HashSet<String> expected=new HashSet<String>();
		// peptides from standard sequence
		expected.add("LLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("GSPAINVAVHVFR");
		expected.add("AADDTWEPFASGK");
		expected.add("TSESGELHGLTTEEEFVEGIYK");
		expected.add("ALGISPFHEHAEVVFTANDSGPR");
		expected.add("YTIAALLSPYSYSTTAVVTNPK");
		// peptide from standard sequence with 1 miss cleavage

		expected.add("MASHRLLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("ASHRLLLLCLAGLVFVSEAGPTGTGESK"); // handle N-term M
		expected.add("LLLLCLAGLVFVSEAGPTGTGESKCPLMVK");
		expected.add("CPLMVKVLDAVR");
		expected.add("VLDAVRGSPAINVAVHVFR");
		expected.add("GSPAINVAVHVFRK");
		expected.add("KAADDTWEPFASGK");
		expected.add("AADDTWEPFASGKTSESGELHGLTTEEEFVEGIYK");
		expected.add("TSESGELHGLTTEEEFVEGIYKVEIDTK");
		expected.add("VEIDTKSYWK");
		expected.add("SYWKALGISPFHEHAEVVFTANDSGPR");
		expected.add("ALGISPFHEHAEVVFTANDSGPRR");
		expected.add("RYTIAALLSPYSYSTTAVVTNPK");
		expected.add("YTIAALLSPYSYSTTAVVTNPKE");
		expected.add("YTIAALLSPYSYSTTAVVTNPKE");

		// simple variant
		expected.add("MASHCLLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("MASHHLLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("ASHCLLLLCLAGLVFVSEAGPTGTGESK");// handle N-term M
		expected.add("ASHHLLLLCLAGLVFVSEAGPTGTGESK");// handle N-term M
		expected.add("LLHLCLAGLVFVSEAGPTGTGESK");
		expected.add("YTIAALLSPYSYTTTAVVTNPK");
		expected.add("YTIAALLSPYSYSTTAVVTNPR");
		expected.add("ALGISPFHEHA");
		expected.add("KPADDTWEPFASGK");
		// with 1 miss cleavage

		expected.add("MGSHRLLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("MASHCLLLLCLAGLVFVSEAGPTGTGESKCPLMVK");
		expected.add("MASHHLLLLCLAGLVFVSEAGPTGTGESKCPLMVK");
		expected.add("MASHRLLHLCLAGLVFVSEAGPTGTGESK");

		expected.add("GSHRLLLLCLAGLVFVSEAGPTGTGESK");// handle N-term M
		expected.add("ASHCLLLLCLAGLVFVSEAGPTGTGESKCPLMVK");// handle N-term M
		expected.add("ASHHLLLLCLAGLVFVSEAGPTGTGESKCPLMVK");// handle N-term M
		expected.add("ASHRLLHLCLAGLVFVSEAGPTGTGESK");// handle N-term M

		expected.add("LLHLCLAGLVFVSEAGPTGTGESKCPLMVK");
		expected.add("SYWKALGISPFHEHA");
		expected.add("RYTIAALLSPYSYTTTAVVTNPK");
		expected.add("YTIAALLSPYSYTTTAVVTNPKE");
		expected.add("RYTIAALLSPYSYSTTAVVTNPR");
		expected.add("YTIAALLSPYSYSTTAVVTNPKD");
		expected.add("KPADDTWEPFASGK");
		expected.add("KPADDTWEPFASGKTSESGELHGLTTEEEFVEGIYK");
		expected.add("YTIAALLSPYSYSTTAVVTNPRE");
		expected.add("YTIAALLSPYSYSTTAVVTNPKD");
		expected.add("GSPAINVAVHVFRKPADDTWEPFASGK");

		// complex variant
		expected.add("THGLTTEEEFVEGIYK");
		expected.add("ALGISPFHEDK");
		expected.add("AEVVFTANDSGPR");
		expected.add("ALGISPFHEAPTAEVVFTANDSGPR");
		expected.add("YTIAALLSPYSYSTTAVTNPK");
		expected.add("LLLLCLAGLVFVSER");
		expected.add("YTIAALLSPY");
		
		
		// with 1 miss cleavage
		expected.add("AADDTWEPFASGKTHGLTTEEEFVEGIYK");
		expected.add("THGLTTEEEFVEGIYKVEIDTK");
		expected.add("ALGISPFHEDKAEVVFTANDSGPR");
		expected.add("SYWKALGISPFHEDK");
		expected.add("AEVVFTANDSGPRR");
		expected.add("ALGISPFHEDKAEVVFTANDSGPR");
		expected.add("SYWKALGISPFHEAPTAEVVFTANDSGPR");
		expected.add("RYTIAALLSPYSYSTTAVTNPK");
		expected.add("YTIAALLSPYSYSTTAVTNPKE");
		expected.add("ALGISPFHEAPTAEVVFTANDSGPRR");
		expected.add("MASHRLLLLCLAGLVFVSER");
		expected.add("ASHRLLLLCLAGLVFVSER");// handle N-term M
		expected.add("LLLLCLAGLVFVSERGSPAINVAVHVFR");
		expected.add("RYTIAALLSPY");
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(fakeTTR.getBytes(StandardCharsets.UTF_8)))), "", "", true, parameters).get(0);
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 40, 1, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry pep : sequences) {
			if (!expected.contains(pep.getSequence())) {
				System.out.println("generated but not expected: "+pep.getSequence());
			}
			expected.remove(pep.getSequence());
		}
		for (String pep : expected) {
			System.out.println("expected but not generated: "+pep);
		}
		assertEquals(0, expected.size());

	}

	public void testMissedCleavageWithVariantsNoNTermMHandling() {
		String fakeTTR=">nxp:NX_P02766-1 \\DbUniqueId=NX_P02766-1 \\PName=Transthyretin isoform Iso 1 \\GName=TTR \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=147 \\SV=266 \\EV=656 \\PE=1 \\ModResPsi=(62|MOD:00041|L-gamma-carboxyglutamic acid)(69|MOD:00047|O-phospho-L-threonine)(72|MOD:00046|O-phospho-L-serine) \\ModRes=(118||N-linked (GlcNAc...)) "
				+"\\VariantSimple=(2|G)(5|H)(5|C)(8|H)(112|*)(137|T)(146|R)(147|D)(147|*)(56|P) "
				+"\\VariantComplex=(70|75|)(110|110|DK)(110|110|APT)(142|142|)(20|40|)(135|147|) \\Processed=(1|20|signal peptide)(21|147|mature protein)"
				+"\nMASHRLLLLCLAGLVFVSEAGPTGTGESKCPLMVKVLDAVRGSPAINVAVHVFRKAADDTWEPFASGKTSESGELHGLTTEEEFVEGIYKVEIDTKSYWKALGISPFHEHAEVVFTANDSGPRRYTIAALLSPYSYSTTAVVTNPKE";
		HashSet<String> expected=new HashSet<String>();
		// peptides from standard sequence
		expected.add("LLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("GSPAINVAVHVFR");
		expected.add("AADDTWEPFASGK");
		expected.add("TSESGELHGLTTEEEFVEGIYK");
		expected.add("ALGISPFHEHAEVVFTANDSGPR");
		expected.add("YTIAALLSPYSYSTTAVVTNPK");
		// peptide from standard sequence with 1 miss cleavage

		expected.add("MASHRLLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("LLLLCLAGLVFVSEAGPTGTGESKCPLMVK");
		expected.add("CPLMVKVLDAVR");
		expected.add("VLDAVRGSPAINVAVHVFR");
		expected.add("GSPAINVAVHVFRK");
		expected.add("KAADDTWEPFASGK");
		expected.add("AADDTWEPFASGKTSESGELHGLTTEEEFVEGIYK");
		expected.add("TSESGELHGLTTEEEFVEGIYKVEIDTK");
		expected.add("VEIDTKSYWK");
		expected.add("SYWKALGISPFHEHAEVVFTANDSGPR");
		expected.add("ALGISPFHEHAEVVFTANDSGPRR");
		expected.add("RYTIAALLSPYSYSTTAVVTNPK");
		expected.add("YTIAALLSPYSYSTTAVVTNPKE");
		expected.add("YTIAALLSPYSYSTTAVVTNPKE");

		// simple variant
		expected.add("MASHCLLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("MASHHLLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("LLHLCLAGLVFVSEAGPTGTGESK");
		expected.add("YTIAALLSPYSYTTTAVVTNPK");
		expected.add("YTIAALLSPYSYSTTAVVTNPR");
		expected.add("ALGISPFHEHA");
		expected.add("KPADDTWEPFASGK");
		// with 1 miss cleavage

		expected.add("MGSHRLLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("MASHCLLLLCLAGLVFVSEAGPTGTGESKCPLMVK");
		expected.add("MASHHLLLLCLAGLVFVSEAGPTGTGESKCPLMVK");
		expected.add("MASHRLLHLCLAGLVFVSEAGPTGTGESK");
		expected.add("LLHLCLAGLVFVSEAGPTGTGESKCPLMVK");
		expected.add("SYWKALGISPFHEHA");
		expected.add("RYTIAALLSPYSYTTTAVVTNPK");
		expected.add("YTIAALLSPYSYTTTAVVTNPKE");
		expected.add("RYTIAALLSPYSYSTTAVVTNPR");
		expected.add("YTIAALLSPYSYSTTAVVTNPKD");
		expected.add("KPADDTWEPFASGK");
		expected.add("KPADDTWEPFASGKTSESGELHGLTTEEEFVEGIYK");
		expected.add("YTIAALLSPYSYSTTAVVTNPRE");
		expected.add("YTIAALLSPYSYSTTAVVTNPKD");
		expected.add("GSPAINVAVHVFRKPADDTWEPFASGK");

		// complex variant
		expected.add("THGLTTEEEFVEGIYK");
		expected.add("ALGISPFHEDK");
		expected.add("AEVVFTANDSGPR");
		expected.add("ALGISPFHEAPTAEVVFTANDSGPR");
		expected.add("YTIAALLSPYSYSTTAVTNPK");
		expected.add("LLLLCLAGLVFVSER");
		expected.add("YTIAALLSPY");


		// with 1 miss cleavage
		expected.add("AADDTWEPFASGKTHGLTTEEEFVEGIYK");
		expected.add("THGLTTEEEFVEGIYKVEIDTK");
		expected.add("ALGISPFHEDKAEVVFTANDSGPR");
		expected.add("SYWKALGISPFHEDK");
		expected.add("AEVVFTANDSGPRR");
		expected.add("ALGISPFHEDKAEVVFTANDSGPR");
		expected.add("SYWKALGISPFHEAPTAEVVFTANDSGPR");
		expected.add("RYTIAALLSPYSYSTTAVTNPK");
		expected.add("YTIAALLSPYSYSTTAVTNPKE");
		expected.add("ALGISPFHEAPTAEVVFTANDSGPRR");
		expected.add("MASHRLLLLCLAGLVFVSER");
		expected.add("LLLLCLAGLVFVSERGSPAINVAVHVFR");
		expected.add("RYTIAALLSPY");
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(fakeTTR.getBytes(StandardCharsets.UTF_8)))), "", "", true, parameters).get(0);
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 40, 1, false, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
		assertEquals(expected.size(), sequences.size());
		for (FastaPeptideEntry pep : sequences) {
			if (!expected.contains(pep.getSequence())) {
				System.out.println("generated but not expected: "+pep.getSequence());
			}
			expected.remove(pep.getSequence());
		}
		for (String pep : expected) {
			System.out.println("expected but not generated: "+pep);
		}
		assertEquals(0, expected.size());

	}

	public void testExpectedModificationsWithVariantsNoNTermMHandling() {
		String fakeTTR=">nxp:NX_Q96B36-1 \\DbUniqueId=NX_Q96B36-1 \\PName=Proline-rich AKT1 substrate 1 isoform Iso 1 \\GName=AKT1S1 \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=256 \\SV=135 \\EV=428 \\PE=1 "
				+ "\\ModResPsi=(3|MOD:00046|O-phospho-L-serine)(51|MOD:00078|omega-N-methyl-L-arginine)(73|MOD:00047|O-phospho-L-threonine)(90|MOD:00047|O-phospho-L-threonine)(97|MOD:00047|O-phospho-L-threonine)(116|MOD:00046|O-phospho-L-serine)(187|MOD:00046|O-phospho-L-serine)(198|MOD:00047|O-phospho-L-threonine)(247|MOD:00046|O-phospho-L-serine)(88|MOD:00046|O-phospho-L-serine)(92|MOD:00046|O-phospho-L-serine)(183|MOD:00046|O-phospho-L-serine)(202|MOD:00046|O-phospho-L-serine)(203|MOD:00046|O-phospho-L-serine)(211|MOD:00046|O-phospho-L-serine)(212|MOD:00046|O-phospho-L-serine)(246|MOD:00047|O-phospho-L-threonine) "
				//+ "\\VariantSimple=(1|V)(3|L)(5|H)(12|T)(13|M)(19|C)(21|W)(21|Q)(23|Q)(26|M)(27|D)(33|T)(40|H)(41|L)(42|S)(43|S)(46|H)(46|C)(47|P)(51|Q)(51|L)(55|V)(59|H)(60|C)(60|H)(61|Y)(63|R)(68|G)(76|Q)(76|W)(78|A)(79|V)(80|S)(82|Q)(86|L)(87|S)(88|C)(88|R)(93|T)(95|W)(95|Q)(95|P)(97|A)(99|V)(103|D)(106|H)(107|K)(114|G)(115|I)(122|T)(122|V)(129|I)(130|A)(135|T)(136|I)(141|H)(141|S)(142|R)(142|H)(142|S)(145|K)(147|E)(149|K)(151|I)(161|S)(162|T)(163|D)(166|I)(168|L)(170|H)(172|V)(174|T)(182|E)(185|S)(186|L)(188|L)(190|F)(190|I)(192|S)(197|T)(198|I)(201|W)(201|Q)(204|N)(207|D)(209|R)(209|L)(214|Y)(215|P)(223|H)(233|A)(233|S)(234|H)(239|V)(242|L)(245|K)(117|R)(212|L)(249|L)(252|Q)(7|K)(51|*)(53|S)(58|V)(64|N)(104|K)(158|*)(163|A)(199|*)(200|V) \\VariantComplex=(39|39|PP)(206|206|)(114|115|) \\Processed=(1|256|mature protein)"
				+"\n"+"MASGRPEELWEAVVGAAERFRARTGTELVLLTAAPPPPPRPGPCAYAAHGRGALAEAARR" +
				"CLHDIALAHRAATAARPPAPPPAPQPPSPTPSPPRPTLAREDNEEDEDEPTETETSGEQL" +
				"GISDNGGLFVMDEDATLQDLPPFCESDPESTDDGSLSEETPAGPPTCSVPPASALPTQQY" +
				"AKSLPVSVPVWGFKEKRTEARSSDEENGPPSSPDLDRIAASMRALVLREAEDTQVFGDLP" +
				"RPRLNTSDFQKLKRKY";
		HashMap<String, String> paramMap=SearchParameterParser.getDefaultParameters();
		paramMap.put("-variable", "S=79.966331,T=79.966331,Y=79.966331");
		SearchParameters parameters=PecanParameterParser.parseParameters(paramMap);
		FastaEntryInterface entry=FastaReader.readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(fakeTTR.getBytes(StandardCharsets.UTF_8)))), "", "", true, parameters).get(0);
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");

		// WITH FIXED AND VARIABLE MODS
		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 4, 40, 0, false, parameters.getAAConstants(), false);
		int haveSTYButNoMods=0;
		for (FastaPeptideEntry fastaPeptideEntry : sequences) {
			if (fastaPeptideEntry.getSequence().indexOf("[79")==-1) {
				if (fastaPeptideEntry.getSequence().indexOf("S")>=0||fastaPeptideEntry.getSequence().indexOf("T")>=0||fastaPeptideEntry.getSequence().indexOf("Y")>=0) {
					haveSTYButNoMods++;
				}
			}
		}
		assertEquals(22, sequences.size());
		assertEquals(9, haveSTYButNoMods);

		sequences=enzyme.digestProtein(entry, 4, 40, 0, false, parameters.getAAConstants(), true);
		for (FastaPeptideEntry fastaPeptideEntry : sequences) {
			assertTrue(fastaPeptideEntry.getSequence().indexOf("[79")>=0);
		}
		assertEquals(10, sequences.size());
	}

	public void testExpectedModificationsWithVariants() {
		String fakeTTR=">nxp:NX_Q96B36-1 \\DbUniqueId=NX_Q96B36-1 \\PName=Proline-rich AKT1 substrate 1 isoform Iso 1 \\GName=AKT1S1 \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=256 \\SV=135 \\EV=428 \\PE=1 "
						+ "\\ModResPsi=(3|MOD:00046|O-phospho-L-serine)(51|MOD:00078|omega-N-methyl-L-arginine)(73|MOD:00047|O-phospho-L-threonine)(90|MOD:00047|O-phospho-L-threonine)(97|MOD:00047|O-phospho-L-threonine)(116|MOD:00046|O-phospho-L-serine)(187|MOD:00046|O-phospho-L-serine)(198|MOD:00047|O-phospho-L-threonine)(247|MOD:00046|O-phospho-L-serine)(88|MOD:00046|O-phospho-L-serine)(92|MOD:00046|O-phospho-L-serine)(183|MOD:00046|O-phospho-L-serine)(202|MOD:00046|O-phospho-L-serine)(203|MOD:00046|O-phospho-L-serine)(211|MOD:00046|O-phospho-L-serine)(212|MOD:00046|O-phospho-L-serine)(246|MOD:00047|O-phospho-L-threonine) "
						//+ "\\VariantSimple=(1|V)(3|L)(5|H)(12|T)(13|M)(19|C)(21|W)(21|Q)(23|Q)(26|M)(27|D)(33|T)(40|H)(41|L)(42|S)(43|S)(46|H)(46|C)(47|P)(51|Q)(51|L)(55|V)(59|H)(60|C)(60|H)(61|Y)(63|R)(68|G)(76|Q)(76|W)(78|A)(79|V)(80|S)(82|Q)(86|L)(87|S)(88|C)(88|R)(93|T)(95|W)(95|Q)(95|P)(97|A)(99|V)(103|D)(106|H)(107|K)(114|G)(115|I)(122|T)(122|V)(129|I)(130|A)(135|T)(136|I)(141|H)(141|S)(142|R)(142|H)(142|S)(145|K)(147|E)(149|K)(151|I)(161|S)(162|T)(163|D)(166|I)(168|L)(170|H)(172|V)(174|T)(182|E)(185|S)(186|L)(188|L)(190|F)(190|I)(192|S)(197|T)(198|I)(201|W)(201|Q)(204|N)(207|D)(209|R)(209|L)(214|Y)(215|P)(223|H)(233|A)(233|S)(234|H)(239|V)(242|L)(245|K)(117|R)(212|L)(249|L)(252|Q)(7|K)(51|*)(53|S)(58|V)(64|N)(104|K)(158|*)(163|A)(199|*)(200|V) \\VariantComplex=(39|39|PP)(206|206|)(114|115|) \\Processed=(1|256|mature protein)"
						+"\n"+"MASGRPEELWEAVVGAAERFRARTGTELVLLTAAPPPPPRPGPCAYAAHGRGALAEAARR" + 
						"CLHDIALAHRAATAARPPAPPPAPQPPSPTPSPPRPTLAREDNEEDEDEPTETETSGEQL" + 
						"GISDNGGLFVMDEDATLQDLPPFCESDPESTDDGSLSEETPAGPPTCSVPPASALPTQQY" + 
						"AKSLPVSVPVWGFKEKRTEARSSDEENGPPSSPDLDRIAASMRALVLREAEDTQVFGDLP" + 
						"RPRLNTSDFQKLKRKY";
		HashMap<String, String> paramMap=SearchParameterParser.getDefaultParameters();
		paramMap.put("-variable", "S=79.966331,T=79.966331,Y=79.966331");
		SearchParameters parameters=PecanParameterParser.parseParameters(paramMap);
		FastaEntryInterface entry=FastaReader.readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(fakeTTR.getBytes(StandardCharsets.UTF_8)))), "", "", true, parameters).get(0);
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		
		// WITH FIXED AND VARIABLE MODS
		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 4, 40, 0, parameters.getAAConstants(), false);
		int haveSTYButNoMods=0;
		for (FastaPeptideEntry fastaPeptideEntry : sequences) {
			if (fastaPeptideEntry.getSequence().indexOf("[79")==-1) {
				if (fastaPeptideEntry.getSequence().indexOf("S")>=0||fastaPeptideEntry.getSequence().indexOf("T")>=0||fastaPeptideEntry.getSequence().indexOf("Y")>=0) {
					haveSTYButNoMods++;
				}
 			}
		}
		assertEquals(24, sequences.size());
		assertEquals(10, haveSTYButNoMods);

		sequences=enzyme.digestProtein(entry, 4, 40, 0, parameters.getAAConstants(), true);
		for (FastaPeptideEntry fastaPeptideEntry : sequences) {
			assertTrue(fastaPeptideEntry.getSequence().indexOf("[79")>=0);
		}
		assertEquals(11, sequences.size());
	}

	public void testModificationsWithVariantsNoNTermMHandling() {
		String fakeTTR=">nxp:NX_P02766-1 \\DbUniqueId=NX_P02766-1 \\PName=Transthyretin isoform Iso 1 \\GName=TTR \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=147 \\SV=266 \\EV=656 \\PE=1 "
				+"\\VariantSimple=(5|C)(8|W)(12|*)(18|P) "
				+"\\VariantComplex=(2|4|)(11|11|WK)(14|14|W*) \\Processed=(1|20|signal peptide)(21|147|mature protein)"
				+"\nMAWHLLWLCLAGLVFVREAGPTGTGESKCPLMV";
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(fakeTTR.getBytes(StandardCharsets.UTF_8)))), "", "", true, parameters).get(0);
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");

		// WITH FIXED MODS
		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 30, 2, false,
				new AminoAcidConstants(new TCharDoubleHashMap(new char[] {'C'}, new double[] {57.0214635}), new ModificationMassMap()), false);

		//123456789012345678901234567890123
		//MAWHLLWLCLAGLVFVREAGPTGTGESKCPLMV

		HashSet<String> expected=getExpectedForModsWithVariantsNoNTermMHandling();

		for (FastaPeptideEntry sequence : sequences) {
			if (!expected.contains(sequence.getSequence())) System.out.println(sequence+" missing from expected");
			assertTrue(expected.remove(sequence.getSequence()));
		}
		assertTrue(expected.size()==0);
	}

	
	public void testModificationsWithVariants() {
		String fakeTTR=">nxp:NX_P02766-1 \\DbUniqueId=NX_P02766-1 \\PName=Transthyretin isoform Iso 1 \\GName=TTR \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=147 \\SV=266 \\EV=656 \\PE=1 "
				+"\\VariantSimple=(5|C)(8|W)(12|*)(18|P) "
				+"\\VariantComplex=(2|4|)(11|11|WK)(14|14|W*) \\Processed=(1|20|signal peptide)(21|147|mature protein)"
				+"\nMAWHLLWLCLAGLVFVREAGPTGTGESKCPLMV";
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FastaEntryInterface entry=FastaReader.readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(fakeTTR.getBytes(StandardCharsets.UTF_8)))), "", "", true, parameters).get(0);
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		
		// WITH FIXED MODS
		ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(entry, 8, 30, 2, new AminoAcidConstants(new TCharDoubleHashMap(new char[] {'C'}, new double[] {57.0214635}), new ModificationMassMap()), false);
		
		//123456789012345678901234567890123
		//MAWHLLWLCLAGLVFVREAGPTGTGESKCPLMV
		
		HashSet<String> expected=getExpectedForModsWithVariants();
		
		for (FastaPeptideEntry sequence : sequences) {
			if (!expected.contains(sequence.getSequence())) System.out.println(sequence+" missing from expected");
			assertTrue(expected.remove(sequence.getSequence()));
		}
		assertTrue(expected.size()==0);
	}

	public HashSet<String> getExpectedForModsWithVariantsNoNTermMHandling() {
		HashSet<String> expected=new HashSet<>();
		//////////////123456789012345678901234567890
		
		expected.add("MAWHLLWLC[57.0214635]LAGLVFVR"); // canonical
		expected.add("MAWHC[57.0214635]LWLC[57.0214635]LAGLVFVR"); //(5|C)
		expected.add("MAWHLLWWC[57.0214635]LAGLVFVR"); //(8|W)
		expected.add("MAWHLLWLC[57.0214635]LA"); //(12|*)
		expected.add("MLLWLC[57.0214635]LAGLVFVR"); //(2|4|)
		expected.add("MAWHLLWLC[57.0214635]LWK"); // (11|11|WK)
		expected.add("MAWHLLWLC[57.0214635]LWKGLVFVR"); // (11|11|WK)
		expected.add("MAWHLLWLC[57.0214635]LAGLW"); // (14|14|W*)

		//////////////123456789012345678901234567890
		expected.add("MAWHLLWLC[57.0214635]LAGLVFVREAGPTGTGESK"); // canonical
		expected.add("MAWHC[57.0214635]LWLC[57.0214635]LAGLVFVREAGPTGTGESK"); //(5|C)
		expected.add("MAWHLLWWC[57.0214635]LAGLVFVREAGPTGTGESK"); //(8|W)
		expected.add("GLVFVREAGPTGTGESK"); //(11|11|WK)
		expected.add("MLLWLC[57.0214635]LAGLVFVREAGPTGTGESK"); //(2|4|)
		//expected.add("VFVREAGPTGTGESK"); // SHOULD NOT EXIST (14|14|W*)
		expected.add("MAWHLLWLC[57.0214635]LAGLVFVRPAGPTGTGESK"); //(18|P)
		
		//////////////8901234567890123456
		expected.add("EAGPTGTGESK"); // canonical
		//expected.add("PAGPTGTGESK"); // SHOULD NOT EXIST (18|P)
		expected.add("EAGPTGTGESKC[57.0214635]PLMV"); // canonical

		// two cleavages
		expected.add("MLLWLC[57.0214635]LAGLVFVREAGPTGTGESKC[57.0214635]PLMV"); // canonical
		expected.add("MAWHLLWLC[57.0214635]LWKGLVFVREAGPTGTGESK"); //(11|11|WK)
		expected.add("GLVFVREAGPTGTGESKC[57.0214635]PLMV"); //(11|11|WK)
		return expected;
	}

	public HashSet<String> getExpectedForModsWithVariants() {
		HashSet<String> expected=new HashSet<>();
		//////////////123456789012345678901234567890

		expected.add("MAWHLLWLC[57.0214635]LAGLVFVR"); // canonical
		expected.add("MAWHC[57.0214635]LWLC[57.0214635]LAGLVFVR"); //(5|C)
		expected.add("MAWHLLWWC[57.0214635]LAGLVFVR"); //(8|W)
		expected.add("MAWHLLWLC[57.0214635]LA"); //(12|*)
		expected.add("MLLWLC[57.0214635]LAGLVFVR"); //(2|4|)
		expected.add("MAWHLLWLC[57.0214635]LWK"); // (11|11|WK)
		expected.add("MAWHLLWLC[57.0214635]LWKGLVFVR"); // (11|11|WK)
		expected.add("MAWHLLWLC[57.0214635]LAGLW"); // (14|14|W*)


		expected.add("AWHLLWLC[57.0214635]LAGLVFVR"); // canonical
		expected.add("AWHC[57.0214635]LWLC[57.0214635]LAGLVFVR"); //(5|C)
		expected.add("AWHLLWWC[57.0214635]LAGLVFVR"); //(8|W)
		expected.add("AWHLLWLC[57.0214635]LA"); //(12|*)
		expected.add("LLWLC[57.0214635]LAGLVFVR"); //(2|4|)
		expected.add("AWHLLWLC[57.0214635]LWK"); // (11|11|WK)
		expected.add("AWHLLWLC[57.0214635]LWKGLVFVR"); // (11|11|WK)
		expected.add("AWHLLWLC[57.0214635]LAGLW"); // (14|14|W*)

		//////////////123456789012345678901234567890
		expected.add("MAWHLLWLC[57.0214635]LAGLVFVREAGPTGTGESK"); // canonical
		expected.add("MAWHC[57.0214635]LWLC[57.0214635]LAGLVFVREAGPTGTGESK"); //(5|C)
		expected.add("MAWHLLWWC[57.0214635]LAGLVFVREAGPTGTGESK"); //(8|W)
		expected.add("GLVFVREAGPTGTGESK"); //(11|11|WK)
		expected.add("MLLWLC[57.0214635]LAGLVFVREAGPTGTGESK"); //(2|4|)
		//expected.add("VFVREAGPTGTGESK"); // SHOULD NOT EXIST (14|14|W*)
		expected.add("MAWHLLWLC[57.0214635]LAGLVFVRPAGPTGTGESK"); //(18|P)

		expected.add("AWHLLWLC[57.0214635]LAGLVFVREAGPTGTGESK"); // canonical // handle N-term M
		expected.add("AWHC[57.0214635]LWLC[57.0214635]LAGLVFVREAGPTGTGESK"); //(5|C) // handle N-term M
		expected.add("AWHLLWWC[57.0214635]LAGLVFVREAGPTGTGESK"); //(8|W) // handle N-term M
		expected.add("LLWLC[57.0214635]LAGLVFVREAGPTGTGESK"); //(2|4|) // handle N-term M
		//expected.add("VFVREAGPTGTGESK"); // SHOULD NOT EXIST (14|14|W*)
		expected.add("AWHLLWLC[57.0214635]LAGLVFVRPAGPTGTGESK"); //(18|P) // handle N-term M

		//////////////8901234567890123456
		expected.add("EAGPTGTGESK"); // canonical
		//expected.add("PAGPTGTGESK"); // SHOULD NOT EXIST (18|P)
		expected.add("EAGPTGTGESKC[57.0214635]PLMV"); // canonical

		// two cleavages
		expected.add("MLLWLC[57.0214635]LAGLVFVREAGPTGTGESKC[57.0214635]PLMV"); // canonical
		expected.add("LLWLC[57.0214635]LAGLVFVREAGPTGTGESKC[57.0214635]PLMV"); // canonical // handle N-term M
		expected.add("MAWHLLWLC[57.0214635]LWKGLVFVREAGPTGTGESK"); //(11|11|WK)
		expected.add("AWHLLWLC[57.0214635]LWKGLVFVREAGPTGTGESK"); //(11|11|WK) // handle N-term M
		expected.add("GLVFVREAGPTGTGESKC[57.0214635]PLMV"); //(11|11|WK)
		return expected;
	}

	public void testXTandemCode() {
		assertEquals("[KR]|{P}", DigestionEnzyme.getEnzyme("Trypsin").toXTandemCode());
		assertEquals("[KR]|{}", DigestionEnzyme.getEnzyme("Trypsin/p").toXTandemCode());
		assertEquals("[K]|{P}", DigestionEnzyme.getEnzyme("Lys-C").toXTandemCode());
		assertEquals("{}|[K]", DigestionEnzyme.getEnzyme("Lys-N").toXTandemCode());
		assertEquals("[R]|{P}", DigestionEnzyme.getEnzyme("Arg-C").toXTandemCode());
		assertEquals("[DE]|{P}", DigestionEnzyme.getEnzyme("Glu-C").toXTandemCode());
		assertEquals("[FWY]|{P}", DigestionEnzyme.getEnzyme("Chymotrypsin").toXTandemCode());
		assertEquals("[FL]|{}", DigestionEnzyme.getEnzyme("Pepsin A").toXTandemCode());
		assertEquals("[AV]|{}", DigestionEnzyme.getEnzyme("Elastase").toXTandemCode());
		assertEquals("{DE}|[AFILMV]", DigestionEnzyme.getEnzyme("Thermolysin").toXTandemCode());
		assertEquals("[]|[]", DigestionEnzyme.getEnzyme("No Enzyme").toXTandemCode());
		assertEquals("{}|{}", new DigestionEnzyme("Non-Specific", "nonspecific", new TCharHashSet(AAs), new TCharHashSet(AAs)).toXTandemCode());
	}

	public void testAvailableEnzymes() {

		final List<DigestionEnzyme> availableEnzymes = getAvailableEnzymes();

		int count=0;

		for (DigestionEnzyme availableEnzyme : availableEnzymes) {
			try {
				final DigestionEnzyme oldEnzyme = oldGetEnzyme(availableEnzyme.getName());
				assertEquals(oldEnzyme, availableEnzyme);
				count++;
			} catch (EncyclopediaException ee) {
				continue;
			}
		}

		assertEquals(14, count);
	}

	// used to test that enzymes didn't change in refactor
	private static DigestionEnzyme oldGetEnzyme(String enzymeName) {
		TCharHashSet n=new TCharHashSet();
		TCharHashSet c=new TCharHashSet();
		if ("Trypsin".equalsIgnoreCase(enzymeName)) {
			n.add('K');
			n.add('R');
			c.addAll(AAs);
			c.remove('P');

			return new DigestionEnzyme("Trypsin", "trypsin", n, c);

		} else if ("Trypsin/p".equalsIgnoreCase(enzymeName)) {
			n.add('K');
			n.add('R');
			c.addAll(AAs);

			return new DigestionEnzyme("Trypsin/p", "trypsinp", n, c);

		} else if ("No Enzyme".equalsIgnoreCase(enzymeName)) {

			return new DigestionEnzyme("No Enzyme", "no_enzyme", n, c);

		} else if ("None".equalsIgnoreCase(enzymeName)) {

			return new DigestionEnzyme("No Enzyme", "no_enzyme", n, c);

		} else if ("Nonspecific".equalsIgnoreCase(enzymeName) || "Nonspecific Enzyme".equalsIgnoreCase(enzymeName)) {

			n.addAll(AAs);
			c.addAll(AAs);
			return new DigestionEnzyme("Nonspecific Enzyme", "nonspecific_enzyme", n, c);

		} else if ("Lys-C".equalsIgnoreCase(enzymeName)) {
			n.add('K');
			c.addAll(AAs);
			c.remove('P');

			return new DigestionEnzyme("Lys-C", "lys-c", n, c);

		} else if ("Lys-N".equalsIgnoreCase(enzymeName)) {
			n.addAll(AAs);
			c.add('K');

			return new DigestionEnzyme("Lys-N", "lys-n", n, c);

		} else if ("Arg-C".equalsIgnoreCase(enzymeName)) {
			n.add('R');
			c.addAll(AAs);
			c.remove('P');

			return new DigestionEnzyme("Arg-C", "arg-c", n, c);

		} else if ("Glu-C".equalsIgnoreCase(enzymeName)) {
			n.add('D'); //Danielle says not to bother
			n.add('E');
			c.addAll(AAs);
			c.remove('P');

			return new DigestionEnzyme("Glu-C", "glu-c", n, c);

		} else if ("Asp-N".equalsIgnoreCase(enzymeName)) {
			n.addAll(AAs);
			c.add('D');
			c.add('E');

			return new DigestionEnzyme("Asp-N", "asp-n", n, c);

		} else if ("Chymotrypsin".equalsIgnoreCase(enzymeName)) {
			n.add('F');
			n.add('Y');
			n.add('W');
			c.addAll(AAs);
			c.remove('P');

			return new DigestionEnzyme("Chymotrypsin", "chymotrypsin", n, c);

		} else if ("Elastase".equalsIgnoreCase(enzymeName)) {
			n.add('A');
			n.add('V');
			c.addAll(AAs);

			return new DigestionEnzyme("Elastase", "elastase", n, c);

		} else if ("Thermolysin".equalsIgnoreCase(enzymeName)) {
			c.add('A');
			c.add('F');
			c.add('I');
			c.add('L');
			c.add('M');
			c.add('V');
			n.addAll(AAs);
			n.remove('D');
			n.remove('E');

			return new DigestionEnzyme("Thermolysin", "thermolysin", n, c);

		} else if ("Pepsin A".equalsIgnoreCase(enzymeName)) {
			n.add('F');
			n.add('L');
			c.addAll(AAs);

			return new DigestionEnzyme("Pepsin A", "pepsin", n, c);
		} else if ("CNBr".equalsIgnoreCase(enzymeName)) {
			n.add('M');
			c.addAll(AAs);

			return new DigestionEnzyme("CNBr", "cnbr", n, c);
		}

		throw new EncyclopediaException("Unknown digestion enzyme ["+enzymeName+"]");
	}


}
