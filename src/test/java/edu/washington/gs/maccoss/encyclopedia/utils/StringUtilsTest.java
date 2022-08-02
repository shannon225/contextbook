package edu.washington.gs.maccoss.encyclopedia.utils;

import junit.framework.TestCase;

public class StringUtilsTest extends TestCase {
	public void testCommon() {
		String[] names=new String[] {"bcs_20161109_yeast_1_to_0_3mz_rep1.mzML", "bcs_20161109_yeast_1_to_0_3mz_rep2.mzML", "bcs_20161109_yeast_1_to_0_3mz_rep3.mzML",
				"bcs_20161109_yeast_1_to_1_3mz_rep1.mzML", "bcs_20161109_yeast_1_to_1_3mz_rep2.mzML", "bcs_20161109_yeast_1_to_1_3mz_rep3.mzML", "bcs_20161109_yeast_1_to_4_3mz_rep1.mzML",
				"bcs_20161109_yeast_1_to_4_3mz_rep2.mzML", "bcs_20161109_yeast_1_to_4_3mz_rep3.mzML", "bcs_20161109_yeast_1_to_9_3mz_rep1.mzML", "bcs_20161109_yeast_1_to_9_3mz_rep2.mzML",
				"bcs_20161109_yeast_1_to_9_3mz_rep3.mzML"};
		
		assertEquals("bcs_20161109_yeast_1_to_.mzML", StringUtils.getCommonName(names, null));
		assertEquals("bcs_20161109_yeast_1_to_combined.mzML", StringUtils.getCommonName(names, "combined"));
	}
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

	public void testTruncationVariable() {
		String[] names=new String[] {"20170426_HZP2N1_01.mzML", "20170426_HZP2N1_02.mzML", "20170426_HZP2N1_03.mzML", "20170426_HZP3N_01.mzML", "20170426_HZP3N_02.mzML", "20170426_HZP3N_03.mzML",
				"20170426_MZP2N1_01.mzML", "20170426_MZP2N1_02.mzML", "20170426_MZP2N1_03.mzML", "20170426_MZP3N_01.mzML", "20170426_MZP3N_02.mzML", "20170426_MZP3N_03.mzML"};
		
		String[] expected=new String[] {"HZP2N1_01", "HZP2N1_02", "HZP2N1_03", "HZP3N_01", "HZP3N_02", "HZP3N_03", "MZP2N1_01", "MZP2N1_02", "MZP2N1_03", "MZP3N_01", "MZP3N_02",
				"MZP3N_03"};
		
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
