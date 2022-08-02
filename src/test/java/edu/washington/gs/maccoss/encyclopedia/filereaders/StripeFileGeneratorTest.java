package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.Iterator;
import java.util.Optional;

import junit.framework.TestCase;

public class StripeFileGeneratorTest extends TestCase {
	public void testFileNameBug() {
		String f="/Volumes/BriansSSD/freezer/121115_bcs_hela_24mz_400_1000_4c_0D.mzML596620318586513071.mzml";

		Optional<String> optional=StripeFileGenerator.getBuggyFileName(f);
		assertTrue(optional.isPresent());
		assertEquals("/Volumes/BriansSSD/freezer/121115_bcs_hela_24mz_400_1000_4c_0D.mzML", optional.get());
	}

	public static void main(String[] args) throws Exception {

		File base=new File("/Users/searleb/Documents/teaching/encyclopedia/final/chromatogram_library");
		File newbase=new File("/Users/searleb/Documents/teaching/encyclopedia/final/test");
		
		File[] origs=new File[] {
				new File(base, "23aug2017_hela_serum_timecourse_wide_1a.dia"),
				new File(base, "23aug2017_hela_serum_timecourse_wide_2f.dia"),
				new File(base, "23aug2017_hela_serum_timecourse_wide_1f.dia"),
				new File(base, "23aug2017_hela_serum_timecourse_wide_3a.dia"),
				new File(base, "23aug2017_hela_serum_timecourse_wide_2a.dia"),
				new File(base, "23aug2017_hela_serum_timecourse_wide_3f.dia")
		};

		File[] news=new File[] {
				new File(newbase, "23aug2017_hela_serum_timecourse_wide_1a.dia"),
				new File(newbase, "23aug2017_hela_serum_timecourse_wide_2f.dia"),
				new File(newbase, "23aug2017_hela_serum_timecourse_wide_1f.dia"),
				new File(newbase, "23aug2017_hela_serum_timecourse_wide_3a.dia"),
				new File(newbase, "23aug2017_hela_serum_timecourse_wide_2a.dia"),
				new File(newbase, "23aug2017_hela_serum_timecourse_wide_3f.dia")
		};

		origs=new File[] {
			new File(base, "23aug2017_hela_serum_timecourse_4mz_narrow_1.dia"),
			new File(base, "23aug2017_hela_serum_timecourse_4mz_narrow_4.dia"),
			new File(base, "23aug2017_hela_serum_timecourse_4mz_narrow_2.dia"),
			new File(base, "23aug2017_hela_serum_timecourse_4mz_narrow_5.dia"),
			new File(base, "23aug2017_hela_serum_timecourse_4mz_narrow_3.dia"),
			new File(base, "23aug2017_hela_serum_timecourse_4mz_narrow_6.dia")
		};
		news=new File[] {
			new File(newbase, "23aug2017_hela_serum_timecourse_4mz_narrow_1.dia"),
			new File(newbase, "23aug2017_hela_serum_timecourse_4mz_narrow_4.dia"),
			new File(newbase, "23aug2017_hela_serum_timecourse_4mz_narrow_2.dia"),
			new File(newbase, "23aug2017_hela_serum_timecourse_4mz_narrow_5.dia"),
			new File(newbase, "23aug2017_hela_serum_timecourse_4mz_narrow_3.dia"),
			new File(newbase, "23aug2017_hela_serum_timecourse_4mz_narrow_6.dia")
		};
		
		for (int i = 0; i < news.length; i++) {
			System.out.println(news[i].getName());
			StripeFile file=(StripeFile)StripeFileGenerator.getFile(origs[i], PecanParameterParser.getDefaultParametersObject());
			file.saveAsFile(news[i]);
			
		}
	}
}
