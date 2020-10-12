package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Optional;

import static org.junit.Assert.assertEquals;

public class MSPReaderTest extends AbstractFileConverterTest {
	private static final String NAME = "MSPReaderTest";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return LibraryFile.DLIB;
	}

	@Test
	public void testSplit() {
		String str="Comment: Single Pep=Tryptic Mods=2/6,C,CAM/16,C,CAM Fullname=R.ALGPAGCEADASAPATCAEMR.C/2 "
				+ "Charge=2 Parent=1053.4561 Se=1(^G1:sc=2.16694e-010) Mz_diff=0ppm Purity=100.0 HCD=44eV "
				+ "Scan=5448 Origfile=\"UNC_P5_whim_8.raw.FT.hcd.ch.MGF\" "
				+ "Nrep=1/2 Sample=\"unc_shared\\Velos_first_runs_UNC_compref\\P5\" "
				+ "Protein=\"gi|54873613|ref|NP_940978.2| agrin precursor\" Unassign_all=0.2243 "
				+ "Unassigned=0.0670 FTResolution=7500 ms2IsolationWidth=2.00 ms1PrecursorAb=387895.41 "
				+ "Precursor1MaxAb=427188.35 PrecursorMonoisoMZ=1053.4556 "
				+ "Filter=\"FTMS + c NSI d Full ms2 1053.96@hcd30.00 [100.00-2000.00]\"";
		
		HashMap<String, String> strings=MSPReader.split(str);

		assertEquals("FTMS + c NSI d Full ms2 1053.96@hcd30.00 [100.00-2000.00]", strings.get("Filter"));
		assertEquals(24, strings.size());
		
		str="Comment: Spec=Consensus Pep=Tryptic Fullname=R.AAAAAAAAAAAAAAAGAGAGAK.Q/2 Mods=0 Parent=798.927 Inst=it Mz_diff=0.003 Mz_exact=798.9268 Mz_av=799.390 Protein=\"sp|P55011|S12A2_HUMAN Solute carrier family 12 member 2 [Homo sapiens]\" Pseq=50/1 Organism=\"human\" Se=4^X62:ex=6.75e-012/2.098e-007,td=7.48e+012/9.807e+013,sd=0/0,hs=90.2/7.974,bs=2.2e-014,b2=6.5e-014,bd=2.27e+015^O80:ex=2.005e-017/1.674e-011,td=3.2e+015/7.262e+017,pr=2.78e-021/2.221e-016,bs=4.16e-021,b2=9.91e-021,bd=1.15e+019^P2:sc=37.3/1.6,dc=24.6/1.6,ps=3.77/0.34,bs=35.7,b2=38.9,bd=23^C9:ex=0.00011/0.0002803,td=455000/8.895e+005,sd=0/0,hs=706/70,bs=1.8e-005,b2=2e-005,bd=2.78e+006 Sample=114/fcrf_anon_00040_none,0,3/fcrf_anon_00041_none,1,4/hlpp_e035_cam,0,1/micro2_01_02_cam,0,1/micro2_04_19_cam,0,1/micro_07_34_cam,0,1/ncta_gastric_mucosa_cam,0,2/pnnl_wqian_hsplasma_none,0,1/resing_hserythroleukq_r_cam,0,1/sri_jurkat_staurosporine_treated_cam,1,2/tcga_a6_3807_01a_22_w_vu_201210_a0218_4i,2,2/tcga_aa_3518_01a_11_w_vu_201209_a0218_3f,3,3/tcga_aa_3526_01a_11_w_vu_201301_a0218_8j,4,4/tcga_aa_3529_01a_12_w_vu_201212_a0218_7b,1,2/tcga_aa_3531_01a_22_w_vu_201301_a0218_8g,6,6/tcga_aa_3534_01a_22_proteome_vu_20130208,0,1/tcga_aa_3552_01a_22_proteome_vu_20130125,1,1/tcga_aa_3554_01a_22_w_vu_201211_a0218_6g,0,1/tcga_aa_3558_01a_22_w_vu_201211_a0218_6a,0,2/tcga_aa_3561_01a_22_w_vu_201208_a0218_1i,1,2/tcga_aa_3664_01a_22_proteome_vu_20130123,2,4/tcga_aa_3666_01a_31_w_vu_201211_a0218_5f,2,3/tcga_aa_3672_01a_22_w_vu_201212_a0218_7j,1,1/tcga_aa_3695_01a_22_proteome_vu_2013013,3,3/tcga_aa_3818_01a_22_proteome_vu_20130210,1,1/tcga_aa_3848_01a_22_w_vu_201212_a0218_7c,4,4/tcga_aa_3864_01a_22_w_vu_201212_a0218_8d,3,3/tcga_aa_3986_01a_12_proteome_vu_20130116,1,1/tcga_aa_3989_01a_22_w_vu_201211_a0218_6b,4,5/tcga_aa_a004_01a_22_w_vu_201210_a0218_4j,0,1/tcga_aa_a00e_01a_31_w_vu_201210_a0218_4e,1,1/tcga_aa_a00f_01a_31_w_vu_201212_a0218_8b,0,1/tcga_aa_a00k_01a_12_w_vu_201212_a0218_7h,1,1/tcga_aa_a00n_01a_32_w_vu_201210_a0218_5d,1,1/tcga_aa_a00n_01a_41_w_vu_201211_a0218_6j,3,3/tcga_aa_a00r_01a_31_w_vu_201212_a0218_8e,1,1/tcga_aa_a00u_01a_41_proteome_vu_20130201,2,2/tcga_aa_a017_01a_22_w_vu_201208_a0218_2b,0,1/tcga_aa_a01d_01a_23_w_vu_201211_a0218_5i,0,1/tcga_aa_a01i_01a_12_w_vu_201211_a0218_6e,5,6/tcga_aa_a01k_01a_31_w_vu_201211_a0218_6d,2,2/tcga_aa_a01p_01a_23_w_vu_201212_a0218_7i,0,1/tcga_aa_a01r_01a_23_w_vu_201208_a0218_2j,1,1/tcga_aa_a01s_01a_23_w_vu_201211_a0218_6f,3,3/tcga_aa_a01t_01a_23_w_vu_201301_a0218_8h,1,3/tcga_aa_a01v_01a_24_w_vu_201211_a0218_6h,5,9/tcga_aa_a01x_01a_23_w_vu_201208_a0218_1h,1,2/tcga_aa_a01z_01a_13_w_vu_201209_a0218_3a,0,3/tcga_aa_a024_01a_32_proteome_vu_20130114,0,2/tcga_aa_a02e_01a_23_w_vu_201209_a0218_3b,0,2/tcga_aa_a02h_01a_32_w_vu_201210_a0218_5b,0,2/tcga_aa_a02j_01a_23_w_vu_201211_a0218_6i,0,2/tcga_aa_a02o_01a_23_proteome_vu_20130205,0,5/tcga_aa_a02r_01a_23_w_vu_201211_a0218_5h,0,1/tcga_aa_a02y_01a_31_w_vu_201209_a0218_3c,0,2/tcga_aa_a03f_01a_41_w_vu_201209_a0218_3i,0,4/tcga_aa_a03j_01a_23_w_vu_201208_a0218_1j,0,3/tcga_af_2692_01a_41_w_vu_201209_a0218_3d,0,2/tcga_af_3913_01a_12_w_vu_201210_a0218_5c,0,1/tcga_ag_3580_01a_22_w_vu_201212_a0218_8c,0,4/tcga_ag_3593_01a_22_w_vu_201209_a0218_3g,0,4/tcga_ag_3594_01a_12_w_vu_201210_a0218_4d,0,2/tcga_ag_a002_01a_23_w_vu_201212_a0218_7g,0,1/tcga_ag_a008_01a_23_proteome_vu_20130127,0,4/tcga_ag_a00c_01a_23_w_vu_201207_a0218_1c,0,3/tcga_ag_a00h_01a_22_w_vu_201212_a0218_7a,0,2/tcga_ag_a00h_01a_31_w_vu_201210_a0218_4b,0,4/tcga_ag_a011_01a_32_w_vu_201208_a0218_2i,0,2/tcga_ag_a015_01a_51_w_vu_201209_a0218_3h,1,6/tcga_ag_a016_01a_23_w_vu_201207_a0218_1d,0,1/tcga_ag_a01l_01a_22_proteome_vu_20130129,0,2/tcga_ag_a01n_01a_23_w_vu_201301_a0218_8i,1,3/tcga_ag_a01w_01a_23_w_vu_201211_a0218_6c,0,3/tcga_ag_a026_01a_71_w_vu_201210_a0218_5a,0,1/tcga_ag_a02n_01a_31_w_vu_201208_a0218_2c,0,2/tcga_ag_a02x_01a_32_proteome_vu_20130213,1,2/tcga_ag_a032_01a_31_w_vu_201211_a0218_5j,0,2/tcga_ag_a036_01a_22_w_vu_201208_a0218_1g,0,1/tum_pxd000125_adenoma_cam,3,3/tum_pxd000125_colon_tissue_cam,1,1/tum_pxd000125_mammary_carcinoma_cam,7,7/tum_pxd000125_stomach_tissue_cam,7,9/uchc_jurkat_depletion_none,0,1/ucsd_06_h293cocl2_total_try_2ul_std_a_200ug_2d34_ltq2_cam,0,2/ucsd_h293_total_try_a_200ug_2d34_081905_ltq1_cam,0,2/ucsd_h293b_total_try_2nd_digest_e_200ug_2d34_new_122305_ltq2_cam,0,1/ucsd_h293b_total_try_2nd_digest_g_200ug_2d34_123005_ltq1_cam,0,1/ucsd_h293b_total_try_b_200ug_2d34_091305_ltq1_cam,0,1/ucsd_h293cocl2_b2_alice_total_try_200ug_a_2d34_120805_ltq2_cam,0,1/ucsd_h293cocl2_total_try_a_200ug_2d34_100605_ltq1_cam,0,3/ucsd_h293cocl2_total_try_b_200ug_2d34_100605_ltq2_cam,0,1/ucsd_hek293_1stbatch_total_trypsin_c_200ug_2umbeads_f17tip_2d34_111105_ltq2_cam,0,1/upenn_postsynaptic_cam,1,2/ut_hek_mudpitvsofgel_cam,0,2/ut_lung_cancer_xenograft_cam,0,4/ut_placenta_hs_villous_cam,0,2/ut_prostate_cancer_t3_cam,0,2/ut_prostate_cancer_t3m_cam,0,3/ut_prostate_cancer_t7_cam,1,4/uta_hek293t_cam,0,2/vu_anon_00169_cam,0,1/vu_anon_00152_cam,0,3/vu_gastric_parietal_cam,0,3/vu_anon_00005_cam,0,2/vu_anon_00035_cam,0,1/vu_anon_00194_cam,0,1/vu_shared_20120605_compref_brplc_update_a2,1,1/vu_shared_20120605_compref_brplc_update_b1,2,2/vu_shared_20120605_compref_brplc_update_b2,3,3/vu_shared_20121219_interstitial_comprefs_whim16_20121013,1,1/vu_anon_00177_cam,0,2/vu_anon_00021_cam,0,1/vu_anon_00024_cam,1,2/vumca_colorectal_cancer_cam,1,4 Nreps=100/266 Missing=0.1792/0.0578 Parent_med=798.92/0.01 Max2med_orig=117.3/41.5 Dotfull=0.780/0.041 Dot_cons=0.854/0.049 Dotbest=0.92 Flags=0,1,82 Unassign_all=0.146 Unassigned=0.021 Naa=22 DUScorr=2.3/1.9/2.9 Dottheory=0.84 Pfin=5.1e+011 Probcorr=20 Tfratio=1.7e+006 Pfract=0";
		strings=MSPReader.split(str);
		assertEquals(32, strings.size());
		
		assertEquals("sp|P55011|S12A2_HUMAN Solute carrier family 12 member 2 [Homo sapiens]", strings.get("Protein"));
		
	}

	@Test
	public void testFastaReader() throws Exception {
		InputStream is=getClass().getResourceAsStream("/truncated.msp");
		ArrayList<LibraryEntry> entries=MSPReader.readMSP(is, "truncated.msp", true);
		assertEquals(4, entries.size());
	}

	//TODO: test readMSP()

	@Test(expected = NullPointerException.class)
	public void testConvertNull() throws Exception {
		MSPReader.convertMSP(null, getFasta().toFile(), out.toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test(expected = FileNotFoundException.class)
	public void testConvertNonexisting() throws Exception {
		final Path msp = Files.createTempFile(tmpDir, NAME, ".msp");
		Files.delete(msp);

		MSPReader.convertMSP(msp.toFile(), getFasta().toFile(), out.toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test
	public void testConvertEmptyFile() throws Exception {
		final Path msp = Files.createTempFile(tmpDir, NAME, ".msp");

		MSPReader.convertMSP(msp.toFile(), getFasta().toFile(), out.toFile(), SearchParameterParser.getDefaultParametersObject());

		final LibraryFile library = new LibraryFile();
		library.openFile(out.toFile());
		try {
			EncyclopediaTestUtils.assertValidDlib(library);
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	Path getFasta() throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/ecoli-190209-contam_correctNL.fasta", tmpDir, NAME, ".fasta");
	}
}
