package edu.washington.gs.maccoss.encyclopedia;

public class EncyclopediaTest {
	public static void main(String[] args) {
		Encyclopedia.main(new String[] { "-i", "/Users/searleb/Documents/projects/encyclopedia/mzml/121115_BCS_HeLa_24mz_400_1000.mzML", 
				//"-l", "/Users/searleb/Documents/projects/encyclopedia/mzml/cptac2_human_hcd_selected.elib",
				"-l", "/Users/searleb/Documents/projects/encyclopedia/mzml/hela_6mz.elib",
				"-targetWindowCenter", "750", 
				"-deconvoluteOverlappingWindows", "true"});
	}
}
