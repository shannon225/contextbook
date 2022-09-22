package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;

import junit.framework.TestCase;

public class TableConcatenatorTest extends TestCase {
	public static void main(String[] args) throws Exception {
//		File output=new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-concatenated.pin");
//		File[] files=new File[] {new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP1_01.pin"),
//				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP2_01.pin"),
//				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP3_01.pin"),
//				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP4_01.pin"),
//				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP5_01.pin"),
//				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP6_01.pin"),
//				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP7_01.pin"),
//				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP8_01.pin"),
//				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRPFT_01.pin"),
//				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRPwash_01.pin")};
		
		File output=new File("/Users/searleb/Documents/whoi/ProteOMZ/library/combined_ProteOMZ_Sept-2018_Metzyme_ORFs_ID_proteins.fasta.trypsin.z3_nce33.spectronaut");
		File[] files=new File[] {
				new File("/Users/searleb/Documents/whoi/ProteOMZ/library/just_z2_ProteOMZ_Sept-2018_Metzyme_ORFs_ID_proteins.fasta.trypsin.z3_nce33.spectronaut"),
				new File("/Users/searleb/Documents/whoi/ProteOMZ/library/just_z3_ProteOMZ_Sept-2018_Metzyme_ORFs_ID_proteins.fasta.trypsin.z3_nce33.spectronaut"),
				new File("/Users/searleb/Documents/whoi/ProteOMZ/library/just_z4_ProteOMZ_Sept-2018_Metzyme_ORFs_ID_proteins.fasta.trypsin.z3_nce33.spectronaut"),
		};
		
		ArrayList<File> list=new ArrayList<>(Arrays.asList(files));
		TableConcatenator.concatenateTables(list, output);
	}

	public void testConcat() throws Exception {
		File f1=FileConcatenatorTest.writeTempFile("/featurefiles/small_23aug2017_hela_serum_timecourse_wide_1a.mzML.features.txt");
		File f2=FileConcatenatorTest.writeTempFile("/featurefiles/small_23aug2017_hela_serum_timecourse_wide_1f.mzML.features.txt");
		File f3=FileConcatenatorTest.writeTempFile("/featurefiles/small_23aug2017_hela_serum_timecourse_wide_2a.mzML.features.txt");
		File f4=FileConcatenatorTest.writeTempFile("/featurefiles/small_23aug2017_hela_serum_timecourse_wide_2f.mzML.features.txt");
		File f5=FileConcatenatorTest.writeTempFile("/featurefiles/small_23aug2017_hela_serum_timecourse_wide_3a.mzML.features.txt");
		File f6=FileConcatenatorTest.writeTempFile("/featurefiles/small_23aug2017_hela_serum_timecourse_wide_3f.mzML.features.txt");
		
		File[] files=new File[] {f1, f2, f3, f4, f5, f6};
		ArrayList<File> list=new ArrayList<>(Arrays.asList(files));
		
		File output=File.createTempFile("output_", ".test");
		
		TableConcatenator.concatenateTables(list, output);
		
		assertEquals(602944, output.length());
		assertEquals(1795, Files.lines(output.toPath()).count());
		
		output.delete();
		
		TableConcatenator.concatenatePINTables(list, output, "primary");
		
		assertEquals(110198, output.length());
		assertEquals(302, Files.lines(output.toPath()).count());
	}

}
