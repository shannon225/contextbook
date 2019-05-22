package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import gnu.trove.map.hash.TCharObjectHashMap;

public class CodonTable {
	public static void main(String[] args) {
		char[] aas="ARNDCQEGHILKMFPSTWYV".toCharArray();
		ArrayList<int[]> matrix=new ArrayList<>();
		matrix.add(new int[] {0,0,0,1523,0,0,957,2074,0,0,0,0,0,0,2330,2936,10938,0,0,9995});
		matrix.add(new int[] {0,0,0,0,11154,13366,0,3645,11820,310,2603,1984,204,0,2329,2214,902,9051,0,0});
		matrix.add(new int[] {0,2,0,2126,0,0,0,0,971,820,0,2547,0,0,0,6179,777,0,556,0});
		matrix.add(new int[] {788,0,6700,0,0,0,2926,3525,2044,0,0,0,0,0,0,0,0,0,2043,1782});
		matrix.add(new int[] {0,2290,0,0,0,0,0,784,0,0,0,0,0,1191,0,1374,0,679,3041,0});
		matrix.add(new int[] {0,3494,0,0,0,0,1840,0,2663,0,616,1264,0,0,1407,0,0,0,0,0});
		matrix.add(new int[] {1382,0,0,3416,0,2570,0,2925,0,0,0,10490,0,0,0,0,0,0,0,1088});
		matrix.add(new int[] {2317,9136,0,4367,1267,0,3635,0,0,0,0,0,0,0,0,6400,0,506,0,3508});
		matrix.add(new int[] {0,2994,582,624,0,1546,0,0,0,0,615,0,0,0,860,0,0,0,2511,0});
		matrix.add(new int[] {0,190,1025,0,0,0,0,0,0,0,1184,178,2111,1090,2,687,5124,0,0,6755});
		matrix.add(new int[] {0,1961,0,0,0,685,0,0,468,1068,0,0,917,4569,5783,1157,0,336,0,4118});
		matrix.add(new int[] {0,3747,2756,0,0,1261,3604,0,0,284,0,0,348,0,0,0,1212,0,0,0});
		matrix.add(new int[] {0,838,0,0,0,0,0,0,0,3035,1413,678,0,0,0,0,3069,0,0,4081});
		matrix.add(new int[] {0,0,0,0,874,0,0,0,0,520,3612,0,0,0,0,1618,0,0,396,698});
		matrix.add(new int[] {2605,2724,0,0,0,716,0,0,991,0,10208,0,0,0,0,6692,2513,0,0,0});
		matrix.add(new int[] {942,2591,3135,0,2243,0,0,2301,0,923,3854,0,0,2863,2391,0,1890,314,995,0});
		matrix.add(new int[] {4456,980,1292,0,0,0,0,0,0,5442,0,796,5325,0,1272,2409,0,0,0,0});
		matrix.add(new int[] {0,1411,0,0,1015,0,0,356,0,0,271,0,0,0,0,382,0,0,0,0});
		matrix.add(new int[] {0,0,468,573,4607,0,0,0,1895,0,0,0,0,573,0,586,0,0,0,0});
		matrix.add(new int[] {3528,0,0,667,0,0,649,1352,0,7460,3781,0,6833,1359,0,0,0,0,0,0});

		AminoAcidConstants aaConstants=new AminoAcidConstants();
		for (int i = 0; i < aas.length; i++) {
			for (int j = 0; j < aas.length; j++) {
				if (matrix.get(i)[j]>0) {
					System.out.println(Math.abs(aaConstants.getMass(aas[i])-aaConstants.getMass(aas[j]))+" "+matrix.get(i)[j]);
				}
			}
		}
	}
	public static void mainSingleAminoAcidSwitch(String[] args) {
		TCharObjectHashMap<String[]> map = new TCharObjectHashMap<>();
		map.put('A', new String[] { "GCT", "GCC", "GCA", "GCG" });
		map.put('R', new String[] { "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" });
		map.put('N', new String[] { "AAT", "AAC" });
		map.put('D', new String[] { "GAT", "GAC" });
		map.put('C', new String[] { "TGT", "TGC" });
		map.put('Q', new String[] { "CAA", "CAG" });
		map.put('E', new String[] { "GAA", "GAG" });
		map.put('G', new String[] { "GGT", "GGC", "GGA", "GGG" });
		map.put('H', new String[] { "CAT", "CAC" });
		map.put('I', new String[] { "ATT", "ATC", "ATA" });
		map.put('L', new String[] { "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" });
		map.put('K', new String[] { "AAA", "AAG" });
		map.put('M', new String[] { "ATG" });
		map.put('F', new String[] { "TTT", "TTC" });
		map.put('P', new String[] { "CCT", "CCC", "CCA", "CCG" });
		map.put('S', new String[] { "TCT", "TCC", "TCA", "TCG", "AGT", "AGC" });
		map.put('T', new String[] { "ACT", "ACC", "ACA", "ACG" });
		map.put('W', new String[] { "TGG" });
		map.put('Y', new String[] { "TAT", "TAC" });
		map.put('V', new String[] { "GTT", "GTC", "GTA", "GTG" });
		
		AminoAcidConstants aaConstants=new AminoAcidConstants();
		for (int i = 0; i < AminoAcidConstants.AAs.length; i++) {
			for (int j = i+1; j < AminoAcidConstants.AAs.length; j++) {
				if (areBitwiseClose(map.get(AminoAcidConstants.AAs[i]), map.get(AminoAcidConstants.AAs[j]))) {
					System.out.println(AminoAcidConstants.AAs[i]+"/"+AminoAcidConstants.AAs[j]+" "+Math.abs(aaConstants.getMass(AminoAcidConstants.AAs[i])-aaConstants.getMass(AminoAcidConstants.AAs[j])));
				}
			}
		}
	}
	
	public static boolean areBitwiseClose(String[] a, String[] b) {
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < b.length; j++) {
				if (areBitwiseClose(a[i], b[j])) return true;
			}
		}
		return false;
	}
	
	public static boolean areBitwiseClose(String a, String b) {
		int diff=areBitwiseDifferent(a, b, 0)+areBitwiseDifferent(a, b, 1)+areBitwiseDifferent(a, b, 2);
		if (diff>1) {
			return false;
		} else {
			return true;
		}
	}

	private static int areBitwiseDifferent(String a, String b, int index) {
		return a.charAt(index)==b.charAt(index)?0:1;
	}
}
