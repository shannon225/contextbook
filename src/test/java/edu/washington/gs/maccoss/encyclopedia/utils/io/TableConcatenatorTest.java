package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

public class TableConcatenatorTest {
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

}
