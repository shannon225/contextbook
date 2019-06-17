package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

public class TableConcatenatorTest {
	public static void main(String[] args) throws Exception {
		File output=new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-concatenated.pin");
		File[] files=new File[] {new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP1_01.pin"),
				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP2_01.pin"),
				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP3_01.pin"),
				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP4_01.pin"),
				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP5_01.pin"),
				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP6_01.pin"),
				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP7_01.pin"),
				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRP8_01.pin"),
				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRPFT_01.pin"),
				new File("/Users/searleb/Downloads/20190130_LUM1_CPBA_EASY04_059_30_SA_1632-01_bRPwash_01.pin")};
		
		ArrayList<File> list=new ArrayList<>(Arrays.asList(files));
		TableConcatenator.concatenateTables(list, output);
	}

}
