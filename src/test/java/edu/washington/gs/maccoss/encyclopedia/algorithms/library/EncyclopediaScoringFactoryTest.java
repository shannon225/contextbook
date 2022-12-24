package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import edu.washington.gs.maccoss.encyclopedia.Batch;

public class EncyclopediaScoringFactoryTest {
	// regression tests
	public static void main(String[] args) {
		String clibPrositTest="/Users/searleb/Documents/encyclopedia/tests/clib_prosit_driver.encxml";
		String clibPanhumanTest="/Users/searleb/Documents/encyclopedia/tests/clib_panhuman_driver.encxml";
		
		Batch.main(new String[] {"-batch", clibPanhumanTest, EncyclopediaScoringFactory.V1_MODE_ARG});
		Batch.main(new String[] {"-batch", clibPrositTest, EncyclopediaScoringFactory.V1_MODE_ARG});

		//Batch.main(new String[] {"-batch", clibPanhumanTest, EncyclopediaScoringFactory.V2_MODE_ARG});
		//Batch.main(new String[] {"-batch", clibPrositTest, EncyclopediaScoringFactory.V2_MODE_ARG});
	}

}
