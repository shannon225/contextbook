package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;

public class IARPATestCase {
	private static final PecanSearchParameters PARAMETERS = new PecanSearchParameters(new AminoAcidConstants(),
			FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"),
			false, true, false);

	private static final File[] files = new File[] {
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_11_1.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_12_2.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_13_3.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_16_4.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_17_5.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_18_6.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_19_7.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_20_8.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_23_9.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_24_12.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_25_13.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_26_14.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_27_15.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_30_16.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_31_17.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_32_18.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_33_19.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_34_20.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_46_21.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_47_22.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_48_23.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_49_24.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_50_25.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_53_26.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_54_27.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_55_28.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_56_29.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_57_30.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_60_31.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_61_32.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_62_33.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_63_34.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_64_35.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_67_36.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_68_37.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_69_38.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_70_39.dia.elib"),
			new File("/Volumes/searle_ssd/localized_individual_results/XXX_2019_0304_RJ_71_40.dia.elib") };

	public static void main(String[] args) throws Exception {
		XCorDIAOneScoringFactory factory = new XCorDIAOneScoringFactory(PARAMETERS);

		File fastaFile = new File(
				"/Volumes/searle_ssd/localized_individual_results/IARPA_var_plus_fasta_20190325.fasta");
		File referenceFile = new File("/Volumes/searle_ssd/localized_individual_results/library/clib.elib");
		LibraryFile reference = new LibraryFile();
		reference.openFile(referenceFile);
		File globalPercolatorOutputFile = new File(
				"/Volumes/searle_ssd/localized_individual_results/2019_quant_reports_concatenated_results.txt");
		File[] sampleFiles = files;

		Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides = PercolatorReader
				.getPassingPeptidesFromTSV(globalPercolatorOutputFile, 0.5f, PARAMETERS.getAAConstants(), false);

		ArrayList<String> peptideList = new ArrayList<String>(Arrays.asList(new String[] { "LDPEAYGSPCAR",
				"LDPEAYGAPCAR", "YTQTFTLHANPAVTYIYNWAYGFGWAATIILIGCAFFFCCLPNYEDDLLGNAK", "YTQTFTLHANR", "VSEALGQGTR",
				"VSDALGQGTR", "FGQGVHHAAGQAGNEAGR", "FGQGVHHAAAQAGNEAGR", "LVWEDTLDK", "LVWEDTLVK", "NTNIAQK",
				"NTNFAQK", "YFSTTEDYDHEITGLR", "YFSTTEDYNHEITGLR", "AEAEALYQTK", "AEAEALYQIK", "ELQNAHNGVNQASK",
				"ELQSAHNGVNQASK", "ALETVQER", "ALETLQER", "FDAPPEAVAAK", "FDAPLEAVAAK", "HGLVATHMLTVR", "HGLVATHTLTVR",
				"AIGGGLSSVGGGSSTIK", "AIGGALSSVGGGSSTIK", "DILTIDIGR", "DILTIDISR", "SWNGSVEILK", "NWNGSVEILK",
				"GSGLGAGQGSNGASVK", "GSGLGAGQGTNGASVK", "WELLQQMNVGTR", "WELLQQMNVDTR", "EAGQFGHDIHHTAGQAGK",
				"EGGQFGHDIHHTAGQAGK", "YQELQITAGR", "YQELQIMAGR",
				"AGGSYGFGGAGSGFGFGGGAGIGFGLGGGAGLAGGFGGPGFPVCPPGGIQEVTVNQSLLTPLNLQIDPAIQR", "AGGSYGFGGAR",
				"QVPGFGVADALGNR", "QVPGFGAADALGNR", "GQHSSGSGQSPGHGQR", "GQHSSGSGQSPGHDQR", "HSASQDGQDTIR",
				"HSAYQDGQDTIR", "QGSSAGSSSSYGQHGSGSR", "QGSSAGSSSSYGPHGSGSR", "GILIDTSR", "GILVDTSR", "HSQSGQGQSAGPR",
				"HSQSGQGQSAGPSTSR", "SALSGHLETVILGLLK", "SALSGHLETLILGLLK", "ELHPVLK", "EFHPVLK", "YGQHGSGSR",
				"YGQHGSGSCQSSGHGR", "GTDECAIESIAVAATPIPK", "GTDECAIESVAVAATPIPK", "VPEPCQPK", "IPEPCQPK", "DLEAHIDSANK",
				"DLEAHVDSANK", "DVTVLQNTDGNNNEAWAK", "DVTVLQNTDGNNNDAWAK", "AEGPEVDVNLPK", "AEDPEVDVNLPK",
				"LISEVDSDGDGEISFQEFLTAAK", "LISEVDSDGDGEISFQEFLTAAR", "VQYDLQK", "VQCDLQK", "SISVSVAGGALLGR",
				"SISVSVAGGALSGR", "GETVSGGNFHGEYPAK", "GETISGGNFHGEYPAK", "FQIATVTEK", "FQTATVTEK",
				"LISEVDSDGDGEISFQEFLTAAK", "LISEVDGDGDGEISFQEFLTAAK", "IMQVVDEK", "VMQVVDEK", "YATTAYMPSEEINLVVK",
				"YATTAYVPSEEINLVVK", "GVALSNVIHK", "GVALSNVVHK", "TYLISSIPLQGAFNYK", "TYLISSIPLHGAFNYK", "HSASQEGQDTIR",
				"HSVSQEGQDTIR", "QGSSAGSSSSYGQHGSGSR", "QGSSAGSSSSCGQHGSGSR", "EWSTFAVGPGHCLQLHDR",
				"EWSTFAVGPGHCLQLNDR", "VPVDVAYR", "VPVDVACR", "QQSHQESTR", "QQSHQESAR", "EELGHLQNDMTSLENDK",
				"EELGHLQNDLTSLENDK", "CLDLGSIIAEVR", "CLDLGSIIAK", "GR", "GHPAVCQPQGR", "QCADLETAIADAEQR",
				"QCAELETAIADAEQR", "VLNDGSVYTAR", "VLNDGTVYTAR", "EQLR", "EQLQQEQALLEEIER", "PEPCISLEPR", "PEPSISLEPR",
				"AASSQTPTMCTTTVTVK", "AASSQTPTMCTTTVTIK", "HSGIGHGQASSAVR", "HAGIGHGQASSAVR", "QGSGSGQSPGHGQR",
				"QGSGSGQSPGHSQR", "VMDVHDGK", "IMDVHDGK", "DITDTLVAVTISEGAHHLDLR", "DITDSLVAVTISEGAHHLDLR",
				"FPSVSLQEASSFFQR", "FPSVSLQEASSFFR", "SSGGSSSVK", "SSGGSSSVR", "QNLEPLFEQYINNLR", "QNMEPLFEQYINNLR",
				"IIPGGIYNADLNDEWVQR", "IILGGIYNADLNDEWVQR" }));
		Collections.sort(peptideList);
		ArrayList<PercolatorPeptide> selectedPeptides = new ArrayList<>();
		for (PercolatorPeptide percolatorPeptide : passingPeptides.x) {
			if (Collections.binarySearch(peptideList, percolatorPeptide.getPeptideSeq()) >= 0) {
				peptideList.remove(percolatorPeptide.getPeptideSeq());
				selectedPeptides.add(percolatorPeptide);
				System.out.println("Found " + percolatorPeptide.getPeptideModSeq() + " at Q-value: "
						+ percolatorPeptide.getQValue());
			}
		}

		System.out.println("Could not find " + peptideList.size() + " peptide(s)");
		for (String missing : peptideList) {
			System.out.println("Missing: " + missing);
		}

		ArrayList<SearchJobData> jobs = new ArrayList<SearchJobData>();
		for (File file : sampleFiles) {
			String absolutePath = file.getAbsolutePath();
			File dia = new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))); // file names are lose
																							// extensions
			XCorDIAJobData job = new XCorDIAJobData(Optional.empty(), dia, fastaFile, factory);
			jobs.add(job);
		}

		ReferencePeakIntegrator.integrateAllPeptides(
				new File("/Volumes/searle_ssd/localized_individual_results/limited_quant_reports.elib"), reference,
				jobs, selectedPeptides, PARAMETERS, new EmptyProgressIndicator());
	}
}
