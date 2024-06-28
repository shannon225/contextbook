package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import junit.framework.TestCase;

public class EncodedAminoAcidTest extends TestCase {
	public void testEncoding() {
		AminoAcidConstants constants=new AminoAcidConstants();
		System.out.println(toString(EncodedAminoAcid.getAAs("Q[-17.026549]AHLC[+57.021464]VLASNC[+57.021464]DEPMYVK", constants)));
		System.out.println(toString(EncodedAminoAcid.getAAs("[+42.010565]S[+79.966331]GSSSVAAM[+15.994915]K", constants)));
		System.out.println(toString(EncodedAminoAcid.getAAs("[+42.010565]SLLDGLASS[+79.966331]PRAPLQSSK", constants)));
		System.out.println(toString(EncodedAminoAcid.getAAs("NDIK[+42.010565]LAAK[+42.010565]LIHTLDDR", constants)));
		System.out.println(toString(EncodedAminoAcid.getAAs("K[+42.010565]PPK[+42.010565]YER", constants)));
		System.out.println(toString(EncodedAminoAcid.getAAs("TLYDESC[+57.0214635]SK[+114.042927]EIQM[+15.994915]AVLLK", constants)));
		System.out.println(toString(EncodedAminoAcid.getAAs("[+229.162932]NVEELNK[+229.162932]", constants)));
		System.out.println(toString(EncodedAminoAcid.getAAs("NVEELNK", constants)));
	}
	
	public String toString(EncodedAminoAcid[] aas) {
		StringBuilder sb=new StringBuilder();
		for (int i = 0; i < aas.length; i++) {
			sb.append(aas[i].toString());
		}
		return sb.toString();
	}
}
