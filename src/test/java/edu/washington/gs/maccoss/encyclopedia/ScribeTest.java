package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;

public class ScribeTest {

	public static void main(String[] args) {
		File dir=new File("/Users/searleb/Documents/dda_library_search/rj_lumos");
		File lib=new File(dir, "uniprot_human_25apr2019.fasta.trypsin.z1-4_nce33.dlib");
		File fasta=new File(dir, "uniprot_human_25apr2019.fasta");
		File[] fs=new File[] {
				//new File(dir, "FU_2016_0627_17_humanHcdLitms2.mzML"),
				// new File(dir, "FU_2016_0627_19_humanCidLitms2.mzML"), 
				//new File(dir, "FU_2016_0627_18_humanHcdOrbims2.dia"), 
				new File(dir, "FU_2016_0627_20_humanCidOrbims2.mzML")
		};

		//lib=new File("/Users/searleb/Documents/dda_library_search/rj_lumos/CIDer/uniprot_human_25apr2019.fasta.trypsin_CIDch2_predictions.dlib");
		//lib=new File("/Users/searleb/Documents/dda_library_search/rj_lumos/CIDer/NIST_CID_whole_human_consensus_final_true_lib.dlib");
		lib=new File("/Users/searleb/Documents/dda_library_search/rj_lumos/CIDer/NIST_CIDer_all-spectra_DBW201129.dlib");
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
				"-ptol", "50", "-ftol", "10", "-lftol", "50"
				//"-ptol", "50", "-ptolunits", "PPM", "-ftol", "1", "-ftolunits", "AMU"
			};
			Scribe.main(scribeArgs);
			duration[i]=System.currentTimeMillis()-startTime;
		}
		
		for (int i = 0; i < duration.length; i++) {
			System.out.println(fs[i].getName()+" \ttotal seconds: "+i/1000f);
		}
	}
}
