package edu.washington.gs.maccoss.encyclopedia;

public class EncyclopediaTest {
	public static void main(String[] args) {
		Encyclopedia.main(new String[] { "-i", "/Users/searleb/Documents/projects/encyclopedia/mzml/yeast/Q_2014_0523_8_0_amol_uL_10mz_overlap.mzML", "-l",
				"/Users/searleb/Documents/projects/encyclopedia/mzml/yeast/yeast_qtof_consensus_final_true_lib.elib",
				//"-targetWindowCenter", "750", 
				"-deconvoluteOverlappingWindows", "true"});
	}
}
