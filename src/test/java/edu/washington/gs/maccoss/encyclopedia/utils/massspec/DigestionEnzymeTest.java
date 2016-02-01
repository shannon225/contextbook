package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import junit.framework.TestCase;

public class DigestionEnzymeTest extends TestCase {
	public void testXanderCleavages() {
		System.out.println("BOLA2T:");
		String sequence="MELSAEYLREKLQRDLEAEHVEVEDTTLNRCSCSFRVLVVSAKFEGKPLLQRHRFCTE";
		for (String name : DigestionEnzyme.getAvailableEnzymes()) {
			DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme(name);
			ArrayList<String> sequences=enzyme.digestProtein(sequence, 6, 40, 0);
			for (String string : sequences) {
				if (string.indexOf("RHRF")>=0) {
					System.out.println(enzyme.getName()+": "+string);
				}
			}
		}
		
		System.out.println("\nBOLA2F:");
		sequence="MELSAEYLREKLQRDLEAEHVEVEDTTLNRCSCSFRVLVVSAKFEGKPLLQRHSLDPSMTIHCDMVITYGLDQLENCQTCGTDYIISVLNLLTLIVEQINTKLPSSFVEKLFIPSSKLLFLRYHKDKEVVAVAHAVYQAMLSLKNIPVLETAYKLILGEMTCALNNLLHSLQLPEACSEIKHEAFKNHVFNVDNAKFVVKFDLSALTTIGNAKNSSL";
		for (String name : DigestionEnzyme.getAvailableEnzymes()) {
			DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme(name);
			ArrayList<String> sequences=enzyme.digestProtein(sequence, 6, 40, 0);
			for (String string : sequences) {
				if (string.indexOf("RHSL")>=0) {
					System.out.println(enzyme.getName()+": "+string);
				}
			}
		}
	}
	
	public void testMissedCleavages() {

		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\n"+"ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF\n"
				+"FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV\n"+"ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK\n"
				+"ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR\n"+"RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE\n"
				+"QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV\n"+"LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL\n"
				+"SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV\n"+"AASQAALGL";

		FastaEntry entry=FastaReader.readFasta(bsa, "").get(0);
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

		ArrayList<String> sequences=enzyme.digestProtein(sequence, 8, 40, 1);
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

		FastaEntry entry=FastaReader.readFasta(bsa, "").get(0);
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
		
		ArrayList<String> sequences=enzyme.digestProtein(sequence, 8, 40, 0);
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

		sequences=enzyme.digestProtein(sequence, 8, 40, 0);
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


		sequences=enzyme.digestProtein(sequence, 8, 40, 0);
		assertEquals(expected.size(), sequences.size());
		for (String peptide : sequences) {
			assertTrue(expected.contains(peptide));
		}
	}

}
