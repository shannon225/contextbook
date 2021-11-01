package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;

public class ScribeTest {

	public static void main(String[] args) {
		File dir=new File("/Users/searleb/Downloads/");
		File lib=new File(dir, "Galaxy6-[TrinityCPTAC.dlib].dlib");
		File fasta=new File(dir, "Galaxy2-[SwissProt_Human_UniProt_Trinity_Validated_cRAP_06102021.fasta].fasta");
		File[] fs=new File[] {
				new File(dir, "Galaxy4-[FN11_N267T272_180min_10ug_C2_111014.mzML].mzml"),
		};

//		File dir=new File("/Users/searleb/Documents/OSU/cptac_grant/cptac_rodland/");
//		File lib=new File(dir, "cptac3_tmt_selected_passed_best.dlib");
//		File fasta=new File(dir, "refseq_human_hg19_refGene_20170329_protein_add_contaminants.fasta");
//		File[] fs=new File[] {
//				new File(dir, "15CPTAC_COprospective_W_PNNL_20170123_B4S3_f01.mzML"),
//				new File(dir, "15CPTAC_COprospective_W_PNNL_20170123_B4S3_f02.mzML"),
//				new File(dir, "15CPTAC_COprospective_W_PNNL_20170123_B4S3_f03.mzML"),
//		};

//		File dir=new File("/Users/searleb/Documents/sperm/");
//		File lib=new File(dir, "psap_albumin_prtc.fasta.trypsin.z3_nce33.dlib");
//		File fasta=new File(dir, "psap_albumin_prtc.fasta");
//		File[] fs=new File[] {
//				new File(dir, "2021mar21_sperm_qc_dda.mzML")
//		};
		
//		File dir=new File("/Users/searleb/Documents/dda_library_search/rj_lumos");
//		File lib=new File(dir, "uniprot_human_25apr2019.fasta.trypsin.z1-4_nce33.dlib");
//		File fasta=new File(dir, "uniprot_human_25apr2019.fasta");
//		File[] fs=new File[] {
//				//new File(dir, "FU_2016_0627_17_humanHcdLitms2.mzML"),
//				// new File(dir, "FU_2016_0627_19_humanCidLitms2.mzML"), 
//				//new File(dir, "FU_2016_0627_18_humanHcdOrbims2.dia"), 
//				new File(dir, "FU_2016_0627_20_humanCidOrbims2.mzML")
//		};

		//lib=new File("/Users/searleb/Documents/dda_library_search/rj_lumos/CIDer/uniprot_human_25apr2019.fasta.trypsin_CIDch2_predictions.dlib");
		//lib=new File("/Users/searleb/Documents/dda_library_search/rj_lumos/CIDer/final_libs/NIST_CID_original__NIST_CID_whole_human_consensus_final_true_lib.dlib");
		//lib=new File("/Users/searleb/Documents/dda_library_search/rj_lumos/CIDer/NIST_CIDer_all-spectra_DBW201129.dlib");
		//lib=new File("/Users/searleb/Documents/dda_library_search/rj_lumos/CIDer/final_libs/NIST_HCD_MS2PIP_PREDICTS__DBW201129.dlib_CIDch2_predictions.dlib");
		//lib=new File("/Users/searleb/Documents/dda_library_search/rj_lumos/CIDer/final_libs/NIST_HCD_CIDer_V2__human_hcd_combined_CIDer_NCE33.dlib");
		//lib=new File("/Users/searleb/Documents/dda_library_search/proteome_tools/ITMS_CID_35_annotated_2019-11-13.dlib");
		//lib=new File("/Users/searleb/Documents/dda_library_search/hela/pan_human_library.dlib");
		//lib=new File(dir, "uniprot_human_25apr2019.fasta.trypsin.z1-4_nce33.dlib");
		//lib=new File("/Users/searleb/Documents/dda_library_search/hela/22oct2017_hela_serum_timecourse_narrow_library.elib");
		//fs=new File[] {new File("/Users/searleb/Documents/dda_library_search/hela/23aug2017_hela_serum_timecourse_pool_dda_001.dia")};
		//lib=new File("/Users/searleb/Documents/dda_library_search/rj_lumos/trypsin_lib.elib");
		
		long[] duration=new long[fs.length];
		for (int i = 0; i < fs.length; i++) {
			long startTime=System.currentTimeMillis();
			String[] scribeArgs=new String[] {
				"-l", lib.getAbsolutePath(), "-i", fs[i].getAbsolutePath(), "-f", fasta.getAbsolutePath(), 
				//"-ptol", "5000", "-ftol", "10", "-lftol", "50", "-percolatorThreshold", "0.05"
				"-ptol", "50", "-ftol", "10", "-lftol", "50"
				//"-ptol", "50", "-ptolunits", "PPM", "-ftol", "10", "-ftolunits", "PPM", "-lftol", "0.6", "-lftolunits", "AMU"
			};
			Scribe.main(scribeArgs);
			duration[i]=System.currentTimeMillis()-startTime;
		}
		
		for (int i = 0; i < duration.length; i++) {
			System.out.println(fs[i].getName()+" \ttotal seconds: "+duration[i]/1000f);
		}
	}
}
