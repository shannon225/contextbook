package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.util.Optional;

import junit.framework.TestCase;

public class StripeFileGeneratorTest extends TestCase {
	public void testFileNameBug() {
		String f="/Volumes/BriansSSD/freezer/121115_bcs_hela_24mz_400_1000_4c_0D.mzML596620318586513071.mzml";

		Optional<String> optional=StripeFileGenerator.getBuggyFileName(f);
		assertTrue(optional.isPresent());
		assertEquals("/Volumes/BriansSSD/freezer/121115_bcs_hela_24mz_400_1000_4c_0D.mzML", optional.get());
	}
}
