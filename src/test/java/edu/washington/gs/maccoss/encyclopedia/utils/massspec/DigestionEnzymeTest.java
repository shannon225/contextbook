package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.ExtendedFastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import gnu.trove.map.hash.TCharDoubleHashMap;
import junit.framework.TestCase;

public class DigestionEnzymeTest extends TestCase {
	public void testReverse() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPENECFLQH\n";

		FastaEntryInterface entry=FastaReader.readFasta(bsa, "").get(0);
		String sequence=entry.getSequence();
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		String reversed=enzyme.reverseProtein(sequence);
		assertEquals("SYASSFLFLLSIFTVWMRFVGRRHADKHAVESRFKFNEEGLDKVHDEFPCQQLYQAFAILVLAKAFETVENVLKDCNEASEDAVCTKDGFLTHLSKLTAVTCLRACCDAMEGYTEKHQLFCENEPEQ", reversed);
	}
	
	public void testFixedMods() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPENECFLQH\n";

		FastaEntryInterface entry=FastaReader.readFasta(bsa, "").get(0);
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
		
		ArrayList<String> sequences=enzyme.digestProtein(entry, 8, 40, 0, constants);
		assertEquals(expected.size(), sequences.size());
		for (String peptide : sequences) {
			assertTrue(expected.contains(peptide));
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
		
		sequences=enzyme.digestProtein(entry, 8, 40, 0, constants);
		assertEquals(expected.size(), sequences.size());
		for (String peptide : sequences) {
			assertTrue(expected.contains(peptide));
		}
	}
	
	public void testModifications() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPENECFLQH\n";

		FastaEntryInterface entry=FastaReader.readFasta(bsa, "").get(0);
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
		
		ArrayList<String> sequences=enzyme.digestProtein(entry, 8, 40, 0, constants);
		assertEquals(expected.size(), sequences.size());
		for (String peptide : sequences) {
			assertTrue(expected.contains(peptide));
		}
	}
	
	public void testXanderCleavages() {
		System.out.println("BOLA2T:");
		String sequence="MELSAEYLREKLQRDLEAEHVEVEDTTLNRCSCSFRVLVVSAKFEGKPLLQRHRFCTE";
		for (DigestionEnzyme enzyme : DigestionEnzyme.getAvailableEnzymes()) {
			ArrayList<String> sequences=enzyme.digestProtein(new FastaEntry(sequence), 6, 40, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
			for (String string : sequences) {
				if (string.indexOf("RHRF")>=0) {
					System.out.println(enzyme.getName()+": "+string);
				}
			}
		}
		
		System.out.println("\nBOLA2F:");
		sequence="MELSAEYLREKLQRDLEAEHVEVEDTTLNRCSCSFRVLVVSAKFEGKPLLQRHSLDPSMTIHCDMVITYGLDQLENCQTCGTDYIISVLNLLTLIVEQINTKLPSSFVEKLFIPSSKLLFLRYHKDKEVVAVAHAVYQAMLSLKNIPVLETAYKLILGEMTCALNNLLHSLQLPEACSEIKHEAFKNHVFNVDNAKFVVKFDLSALTTIGNAKNSSL";
		for (DigestionEnzyme enzyme : DigestionEnzyme.getAvailableEnzymes()) {
			ArrayList<String> sequences=enzyme.digestProtein(new FastaEntry(sequence), 6, 40, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
			for (String string : sequences) {
				if (string.indexOf("RHSL")>=0) {
					System.out.println(enzyme.getName()+": "+string);
				}
			}
		}
	}
	public void testNoEnzyme() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\n"+"ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF\n"
				+"FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV\n"+"ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK\n"
				+"ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR\n"+"RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE\n"
				+"QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV\n"+"LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL\n"
				+"SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV\n"+"AASQAALGL";

		FastaEntryInterface entry=FastaReader.readFasta(bsa, "").get(0);
		String sequence=entry.getSequence();
		
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("no enzyme");

		ArrayList<String> sequences=enzyme.digestProtein(entry, 8, 99999999, 1, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
		assertEquals(1, sequences.size());
		assertEquals(entry.getSequence(), sequences.get(0));
	}
	
	public void testMissedCleavages() {

		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\n"+"ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF\n"
				+"FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV\n"+"ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK\n"
				+"ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR\n"+"RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE\n"
				+"QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV\n"+"LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL\n"
				+"SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV\n"+"AASQAALGL";

		FastaEntryInterface entry=FastaReader.readFasta(bsa, "").get(0);
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

		ArrayList<String> sequences=enzyme.digestProtein(entry, 8, 40, 1, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
		assertEquals(expected.size(), sequences.size());
		for (String peptide : sequences) {
			assertTrue(expected.contains(peptide));
		}
		
	}
	public void testEnzymes() {
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\n"+"ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF\n"
				+"FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV\n"+"ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK\n"
				+"ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR\n"+"RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE\n"
				+"QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV\n"+"LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL\n"
				+"SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV\n"+"AASQAALGL";

		FastaEntryInterface entry=FastaReader.readFasta(bsa, "").get(0);
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
		
		ArrayList<String> sequences=enzyme.digestProtein(entry, 8, 40, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
		assertEquals(expected.size(), sequences.size());
		for (String peptide : sequences) {
			assertTrue(expected.contains(peptide));
		}
		expected.clear();
		
		enzyme=DigestionEnzyme.getEnzyme("Chymotrypsin");
		expected.add("RRDAHKSEVAHRF");
		expected.add("AKTCVADESAENCDKSL");
		expected.add("GEMADCCAKQEPERNECF");
		expected.add("QHKDDNPNLPRL");
		expected.add("VRPEVDVMCTAF");
		expected.add("EIARRHPY");
		expected.add("TECCQAADKAACL");
		expected.add("RDEGKASSAKQRL");
		expected.add("SQRFPKAEF");
		expected.add("TKVHTECCHGDL");
		expected.add("ECADDRADL");
		expected.add("ICENQDSISSKL");
		expected.add("KECCEKPL");
		expected.add("EKSHCIAEVENDEMPADLPSL");
		expected.add("VESKDVCKNY");
		expected.add("EKCCAAADPHECY");
		expected.add("TKKVPQVSTPTL");
		expected.add("GKVGSKCCKHPEAKRMPCAEDY");
		expected.add("HEKTPVSDRVTKCCTESL");
		expected.add("SEKERQIKKQTAL");
		expected.add("VKHKPKATKEQL");
		expected.add("VEKCCKADDKETCF");
		expected.add("VAASQAAL");

		sequences=enzyme.digestProtein(entry, 8, 40, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
		assertEquals(expected.size(), sequences.size());
		for (String peptide : sequences) {
			assertTrue(expected.contains(peptide));
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


		sequences=enzyme.digestProtein(entry, 8, 40, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
		assertEquals(expected.size(), sequences.size());
		for (String peptide : sequences) {
			assertTrue(expected.contains(peptide));
		}
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

		FastaEntryInterface entry=FastaReader.readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(fakeTTR.getBytes(StandardCharsets.UTF_8)))), "", "", true).get(0);

		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");

		ArrayList<String> sequences=enzyme.digestProtein(entry, 8, 40, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
		assertEquals(expected.size(), sequences.size());
		for (String pep : sequences) {
			if (!expected.contains(pep)) {
				System.out.println("generated but not expected: "+pep);
			}
		}
		expected.removeAll(sequences);
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
		FastaEntryInterface entry=FastaReader.readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(fakeTTR.getBytes(StandardCharsets.UTF_8)))), "", "", true).get(0);
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		ArrayList<String> sequences=enzyme.digestProtein(entry, 8, 40, 1, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
		assertEquals(expected.size(), sequences.size());
		for (String pep : sequences) {
			if (!expected.contains(pep)) {
				System.out.println("generated but not expected: "+pep);
			}
		}
		expected.removeAll(sequences);
		for (String pep : expected) {
			System.out.println("expected but not generated: "+pep);
		}
		assertEquals(0, expected.size());

	}

}
