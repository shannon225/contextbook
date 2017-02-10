package edu.washington.gs.maccoss.encyclopedia.utils;

import junit.framework.TestCase;

public class StringUtilsTest extends TestCase {
	public void testTruncation() {
		String[] names=new String[] {"bcs_20161109_yeast_1_to_0_3mz_rep1.mzML", "bcs_20161109_yeast_1_to_0_3mz_rep2.mzML", "bcs_20161109_yeast_1_to_0_3mz_rep3.mzML",
				"bcs_20161109_yeast_1_to_1_3mz_rep1.mzML", "bcs_20161109_yeast_1_to_1_3mz_rep2.mzML", "bcs_20161109_yeast_1_to_1_3mz_rep3.mzML", "bcs_20161109_yeast_1_to_4_3mz_rep1.mzML",
				"bcs_20161109_yeast_1_to_4_3mz_rep2.mzML", "bcs_20161109_yeast_1_to_4_3mz_rep3.mzML", "bcs_20161109_yeast_1_to_9_3mz_rep1.mzML", "bcs_20161109_yeast_1_to_9_3mz_rep2.mzML",
				"bcs_20161109_yeast_1_to_9_3mz_rep3.mzML"};
		
		String[] expected=new String[] {"0_3mz_rep1", "0_3mz_rep2", "0_3mz_rep3", "1_3mz_rep1", "1_3mz_rep2", "1_3mz_rep3", "4_3mz_rep1", "4_3mz_rep2", "4_3mz_rep3", "9_3mz_rep1", "9_3mz_rep2",
				"9_3mz_rep3"};
		
		String[] truncated=StringUtils.getUniquePortion(names);
		for (int i=0; i<truncated.length; i++) {
			assertEquals(expected[i], truncated[i]);
		}

		// test to make sure boundary cases are ok
		String[] doubleTruncated=StringUtils.getUniquePortion(truncated);
		for (int i=0; i<doubleTruncated.length; i++) {
			assertEquals(truncated[i], doubleTruncated[i]);
		}
	}
}
