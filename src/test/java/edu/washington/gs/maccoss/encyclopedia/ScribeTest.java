package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;

public class ScribeTest {

	public static void main(String[] args) {
		File dir=new File("/Users/searleb/Documents/dda_library_search/rj_lumos");
		File lib=new File(dir, "uniprot_human_25apr2019.fasta.z2_nce33.dlib");
		File fasta=new File(dir, "uniprot_human_25apr2019.fasta");
		File[] fs=new File[] {
				// new File(dir, "FU_2016_0627_17_humanHcdLitms2.mzML"),
				// new File(dir, "FU_2016_0627_19_humanCidLitms2.mzML"), 
				new File(dir, "FU_2016_0627_18_humanHcdOrbims2.mzML"), 
				//new File(dir, "FU_2016_0627_20_humanCidOrbims2.mzML")
		};
		
		long[] duration=new long[fs.length];
		for (int i = 0; i < fs.length; i++) {
			long startTime=System.currentTimeMillis();
			String[] scribeArgs=new String[] {
				"-l", lib.getAbsolutePath(), "-i", fs[i].getAbsolutePath(), "-f", fasta.getAbsolutePath(), 
				"-ptol", "25", "-ftol", "25"
			};
			Scribe.main(scribeArgs);
			duration[i]=System.currentTimeMillis()-startTime;
		}
		
		for (int i = 0; i < duration.length; i++) {
			System.out.println(fs[i].getName()+" \ttotal seconds: "+i/1000f);
		}
	}
}
