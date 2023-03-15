package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import edu.washington.gs.maccoss.encyclopedia.Batch;

public class EncyclopediaScoringFactoryTest {
	// regression tests
	public static void main(String[] args) {
//		EncyclopediaScoringFactory.USE_LEGACY_SCORING_SYSTEM=true;
//		String clibPrositTest="/Users/searleb/Documents/encyclopedia/tests/clib_prosit_driver.encxml";
//		String clibPanhumanTest="/Users/searleb/Documents/encyclopedia/tests/clib_panhuman_driver.encxml";
//		
//		Batch.main(new String[] {"-batch", clibPanhumanTest});
//		Batch.main(new String[] {"-batch", clibPrositTest});

		EncyclopediaScoringFactory.USE_LEGACY_SCORING_SYSTEM=false;
		String clibV2PrositTest="/Users/searleb/Documents/encyclopedia/tests/clibV2_prosit_driver.encxml";
		String clibV2PanhumanTest="/Users/searleb/Documents/encyclopedia/tests/clibV2_panhuman_driver.encxml";
		
		Batch.main(new String[] {"-batch", clibV2PanhumanTest});
		Batch.main(new String[] {"-batch", clibV2PrositTest});
	}

}
