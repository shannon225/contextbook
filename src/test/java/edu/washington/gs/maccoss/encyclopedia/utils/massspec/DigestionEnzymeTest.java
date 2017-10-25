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
	public void testDigestionWithExtendedFastaEntry() {
		// adding some random complex variant complex cases into TTR to test
		// digestion function
		String fakeTTR=">nxp:NX_P02766-1 \\DbUniqueId=NX_P02766-1 \\PName=Transthyretin isoform Iso 1 \\GName=TTR \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=147 \\SV=266 \\EV=656 \\PE=1 \\ModResPsi=(62|MOD:00041|L-gamma-carboxyglutamic acid)(69|MOD:00047|O-phospho-L-threonine)(72|MOD:00046|O-phospho-L-serine) \\ModRes=(118||N-linked (GlcNAc...)) "
				+"\\VariantSimple=(2|G)(5|H)(5|C)(8|H)(9|V)(9|F)(13|E)(13|R)(18|F)(18|A)(19|D)(23|M)(26|D)(26|S)(30|R)(32|P)(33|T)(33|I)(38|E)(38|G)(40|I)(41|Q)(42|D)(43|N)(44|S)(45|T)(46|V)(47|S)(48|M)(50|L)(50|G)(50|A)(50|M)(51|R)(53|C)(53|V)(53|L)(53|I)(54|T)(55|N)(55|M)(56|P)(58|A)(58|V)(61|L)(62|D)(62|G)(64|L)(64|S)(65|S)(65|D)(65|T)(67|E)(67|V)(67|R)(67|A)(69|I)(69|A)(70|N)(70|G)(70|R)(70|I)(72|P)(72|Y)(73|E)(74|G)(74|K)(74|S)(75|Q)(75|P)(76|Y)(78|H)(78|R)(79|I)(79|A)(79|K)(80|A)(81|G)(81|K)(82|D)(84|L)(88|L)(89|H)(90|N)(91|A)(93|V)(94|H)(95|S)(97|Y)(98|F)(99|*)(100|E)(104|N)(104|T)(104|L)(104|S)(105|P)(109|K)(109|D)(109|Q)(110|N)(111|S)(112|*)(114|A)(117|G)(117|S)(119|A)(121|D)(121|S)(121|A)(122|R)(123|H)(124|C)(124|H)(125|H)(125|*)(126|N)(126|I)(127|M)(127|V)(129|T)(131|M)(134|H)(134|C)(136|V)(136|H)(136|S)(139|M)(140|V)(140|S)(142|A)(142|I)(144|S)(145|S)(52|E)(80|I)(84|S)(89|C)(90|Q)(101|S)(101|T)(103|A)(118|K)(119|N)(120|P)(123|S)(123|C)(127|T)(129|V)(133|S)(135|F)(136|C)(137|T)(146|R)(147|D)(147|*) "
				+"\\VariantComplex=(3|3|SHAD)(4|4|)(32|32|LM)(70|75|)(110|110|DK)(110|110|APT)(142|142|) \\Processed=(1|20|signal peptide)(21|147|mature protein)\n"
				+"MASHRLLLLCLAGLVFVSEAGPTGTGESKCPLMVKVLDAVRGSPAINVAVHVFRKAADDTWEPFASGKTSESGELHGLTTEEEFVEGIYKVEIDTKSYWKALGISPFHEHAEVVFTANDSGPRRYTIAALLSPYSYSTTAVVTNPKE";

		HashSet<String> expected=new HashSet<String>();
		// peptides from standard sequence
		expected.add("LLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("GSPAINVAVHVFR");
		expected.add("AADDTWEPFASGK");
		expected.add("TSESGELHGLTTEEEFVEGIYK");
		expected.add("ALGISPFHEHAEVVFTANDSGPR");
		expected.add("YTIAALLSPYSYSTTAVVTNPK");

		// additional peptides after considering simple variant annotation
		expected.add("MASHHLLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("MASHCLLLLCLAGLVFVSEAGPTGTGESK");
		expected.add("LLHLCLAGLVFVSEAGPTGTGESK");
		expected.add("LLLVCLAGLVFVSEAGPTGTGESK");
		expected.add("LLLFCLAGLVFVSEAGPTGTGESK");
		expected.add("LLLLCLAELVFVSEAGPTGTGESK");
		expected.add("LLLLCLAR");
		expected.add("LVFVSEAGPTGTGESK");
		expected.add("LLLLCLAGLVFVFEAGPTGTGESK");
		expected.add("LLLLCLAGLVFVAEAGPTGTGESK");
		expected.add("LLLLCLAGLVFVSDAGPTGTGESK");
		expected.add("LLLLCLAGLVFVSEAGPMGTGESK");
		expected.add("LLLLCLAGLVFVSEAGPTGTDESK");
		expected.add("LLLLCLAGLVFVSEAGPTGTSESK");
		expected.add("VLDAVQGSPAINVAVHVFR");
		expected.add("DSPAINVAVHVFR");
		expected.add("GNPAINVAVHVFR");
		expected.add("GSSAINVAVHVFR");
		expected.add("GSPTINVAVHVFR");
		expected.add("GSPAVNVAVHVFR");
		expected.add("GSPAISVAVHVFR");
		expected.add("GSPAINMAVHVFR");
		expected.add("GSPAINVALHVFR");
		expected.add("GSPAINVAGHVFR");
		expected.add("GSPAINVAAHVFR");
		expected.add("GSPAINVAMHVFR");
		expected.add("GSPAINVAVR");
		expected.add("GSPAINVAVHVCR");
		expected.add("GSPAINVAVHVVR");
		expected.add("GSPAINVAVHVLR");
		expected.add("GSPAINVAVHVIR");
		expected.add("GSPAINVAVHVFTK");
		expected.add("NAADDTWEPFASGK");
		expected.add("MAADDTWEPFASGK");
		expected.add("KPADDTWEPFASGK");
		expected.add("AAADTWEPFASGK");
		expected.add("AAVDTWEPFASGK");
		expected.add("AADDTLEPFASGK");
		expected.add("AADDTWDPFASGK");
		expected.add("AADDTWGPFASGK");
		expected.add("AADDTWEPLASGK");
		expected.add("AADDTWEPSASGK");
		expected.add("AADDTWEPFSSGK");
		expected.add("AADDTWEPFDSGK");
		expected.add("AADDTWEPFTSGK");
		expected.add("AADDTWEPFASEK");
		expected.add("AADDTWEPFASVK");
		expected.add("AADDTWEPFASR");
		expected.add("AADDTWEPFASAK");
		expected.add("ISESGELHGLTTEEEFVEGIYK");
		expected.add("ASESGELHGLTTEEEFVEGIYK");
		expected.add("TNESGELHGLTTEEEFVEGIYK");
		expected.add("TGESGELHGLTTEEEFVEGIYK");
		expected.add("ESGELHGLTTEEEFVEGIYK");
		expected.add("TIESGELHGLTTEEEFVEGIYK");
		expected.add("TSEPGELHGLTTEEEFVEGIYK");
		expected.add("TSEYGELHGLTTEEEFVEGIYK");
		expected.add("TSESEELHGLTTEEEFVEGIYK");
		expected.add("TSESGGLHGLTTEEEFVEGIYK");
		expected.add("LHGLTTEEEFVEGIYK");
		expected.add("TSESGSLHGLTTEEEFVEGIYK");
		expected.add("TSESGEQHGLTTEEEFVEGIYK");
		expected.add("TSESGEPHGLTTEEEFVEGIYK");
		expected.add("TSESGELYGLTTEEEFVEGIYK");
		expected.add("TSESGELHGHTTEEEFVEGIYK");
		expected.add("TSESGELHGR");
		expected.add("TTEEEFVEGIYK");
		expected.add("TSESGELHGLITEEEFVEGIYK");
		expected.add("TSESGELHGLATEEEFVEGIYK");
		expected.add("TSESGELHGLK");
		expected.add("TEEEFVEGIYK");
		expected.add("TSESGELHGLTAEEEFVEGIYK");
		expected.add("TSESGELHGLTTGEEFVEGIYK");
		expected.add("TSESGELHGLTTK");
		expected.add("EEFVEGIYK");
		expected.add("TSESGELHGLTTEDEFVEGIYK");
		expected.add("TSESGELHGLTTEEELVEGIYK");
		expected.add("TSESGELHGLTTEEEFVEGLYK");
		expected.add("TSESGELHGLTTEEEFVEGIHK");
		expected.add("TSESGELHGLTTEEEFVEGIYNVEIDTK");
		expected.add("SYWEALGISPFHEHAEVVFTANDSGPR");
		expected.add("ALGNSPFHEHAEVVFTANDSGPR");
		expected.add("ALGTSPFHEHAEVVFTANDSGPR");
		expected.add("ALGLSPFHEHAEVVFTANDSGPR");
		expected.add("ALGSSPFHEHAEVVFTANDSGPR");
		expected.add("ALGIPPFHEHAEVVFTANDSGPR");
		expected.add("ALGISPFHK");
		expected.add("HAEVVFTANDSGPR");
		expected.add("ALGISPFHDHAEVVFTANDSGPR");
		expected.add("ALGISPFHQHAEVVFTANDSGPR");
		expected.add("ALGISPFHENAEVVFTANDSGPR");
		expected.add("ALGISPFHEHSEVVFTANDSGPR");
		expected.add("ALGISPFHEHA");
		expected.add("ALGISPFHEHAEVAFTANDSGPR");
		expected.add("ALGISPFHEHAEVVFTGNDSGPR");
		expected.add("ALGISPFHEHAEVVFTSNDSGPR");
		expected.add("ALGISPFHEHAEVVFTANASGPR");
		expected.add("ALGISPFHEHAEVVFTANDSDPR");
		expected.add("ALGISPFHEHAEVVFTANDSSPR");
		expected.add("ALGISPFHEHAEVVFTANDSAPR");
		expected.add("ALGISPFHEHAEVVFTANDSGR");
		expected.add("ALGISPFHEHAEVVFTANDSGPHR");
		expected.add("YTIAALLSPYSYSTTAVVTNPK");
		expected.add("CYTIAALLSPYSYSTTAVVTNPK");
		expected.add("HYTIAALLSPYSYSTTAVVTNPK");
		expected.add("HTIAALLSPYSYSTTAVVTNPK");
		expected.add("YNIAALLSPYSYSTTAVVTNPK");
		expected.add("YIIAALLSPYSYSTTAVVTNPK");
		expected.add("YTMAALLSPYSYSTTAVVTNPK");
		expected.add("YTVAALLSPYSYSTTAVVTNPK");
		expected.add("YTIATLLSPYSYSTTAVVTNPK");
		expected.add("YTIAALMSPYSYSTTAVVTNPK");
		expected.add("YTIAALLSPHSYSTTAVVTNPK");
		expected.add("YTIAALLSPCSYSTTAVVTNPK");
		expected.add("YTIAALLSPYSVSTTAVVTNPK");
		expected.add("YTIAALLSPYSHSTTAVVTNPK");
		expected.add("YTIAALLSPYSSSTTAVVTNPK");
		expected.add("YTIAALLSPYSYSTMAVVTNPK");
		expected.add("YTIAALLSPYSYSTTVVVTNPK");
		expected.add("YTIAALLSPYSYSTTSVVTNPK");
		expected.add("YTIAALLSPYSYSTTAVATNPK");
		expected.add("YTIAALLSPYSYSTTAVITNPK");
		expected.add("YTIAALLSPYSYSTTAVVTSPK");
		expected.add("YTIAALLSPYSYSTTAVVTNSK");
		expected.add("GSPAINVAVHEFR");
		expected.add("TSESGELHGLTIEEEFVEGIYK");
		expected.add("TSESGELHGLTTEEESVEGIYK");
		expected.add("TSESGELHGLTTEEEFVEGICK");
		expected.add("TSESGELHGLTTEEEFVEGIYQVEIDTK");
		expected.add("SLGISPFHEHAEVVFTANDSGPR");
		expected.add("TLGISPFHEHAEVVFTANDSGPR");
		expected.add("ALAISPFHEHAEVVFTANDSGPR");
		expected.add("ALGISPFHEHAEVVFTAK");
		expected.add("ALGISPFHEHAEVVFTANNSGPR");
		expected.add("ALGISPFHEHAEVVFTANDPGPR");
		expected.add("ALGISPFHEHAEVVFTANDSGPSR");
		expected.add("ALGISPFHEHAEVVFTANDSGPCR");
		expected.add("YTIAALLSPYSYSTTAVVTNPK");
		expected.add("YTTAALLSPYSYSTTAVVTNPK");
		expected.add("YTIAVLLSPYSYSTTAVVTNPK");
		expected.add("YTIAALLSSYSYSTTAVVTNPK");
		expected.add("YTIAALLSPYFYSTTAVVTNPK");
		expected.add("YTIAALLSPYSCSTTAVVTNPK");
		expected.add("YTIAALLSPYSYTTTAVVTNPK");
		expected.add("YTIAALLSPYSYSTTAVVTNPR");

		// additional peptides after considering complex variant annotation
		expected.add("MASHADHR");
		expected.add("THGLTTEEEFVEGIYK");
		expected.add("ALGISPFHEDK");
		expected.add("AEVVFTANDSGPR");
		expected.add("ALGISPFHEAPTAEVVFTANDSGPR");
		expected.add("YTIAALLSPYSYSTTAVTNPK");

		FastaEntryInterface entry=FastaReader.readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(fakeTTR.getBytes(StandardCharsets.UTF_8)))), "", "", true).get(0);
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");

		ArrayList<String> sequences=enzyme.digestProtein(entry, 8, 40, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
		assertEquals(expected.size(), sequences.size());
		expected.removeAll(sequences);
		assertEquals(0, expected.size());
	}
	
	

}
