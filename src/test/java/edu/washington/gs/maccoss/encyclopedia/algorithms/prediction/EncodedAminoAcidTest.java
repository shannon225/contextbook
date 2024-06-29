package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import junit.framework.TestCase;

public class EncodedAminoAcidTest extends TestCase {
	public void testEncoding() {
		AminoAcidConstants constants=new AminoAcidConstants();
		assertEquals("Q[-17],A,H,L,C[+57],V,L,A,S,N,C[+57],D,E,P,M,Y,V,K", toString(EncodedAminoAcid.getAAs("Q[-17.026549]AHLC[+57.021464]VLASNC[+57.021464]DEPMYVK", constants)));
		assertEquals("[+42],S[+80],G,S,S,S,V,A,A,M[+16],K", toString(EncodedAminoAcid.getAAs("[+42.010565]S[+79.966331]GSSSVAAM[+15.994915]K", constants)));
		assertEquals("[+42],S,L,L,D,G,L,A,S,S[+80],P,R,A,P,L,Q,S,S,K", toString(EncodedAminoAcid.getAAs("[+42.010565]SLLDGLASS[+79.966331]PRAPLQSSK", constants)));
		assertEquals("N,D,I,K[+42],L,A,A,K[+42],L,I,H,T,L,D,D,R", toString(EncodedAminoAcid.getAAs("NDIK[+42.010565]LAAK[+42.010565]LIHTLDDR", constants)));
		assertEquals("K[+42],P,P,K[+42],Y,E,R", toString(EncodedAminoAcid.getAAs("K[+42.010565]PPK[+42.010565]YER", constants)));
		assertEquals("[+42],K,P,P,K[+42],Y,E,R", toString(EncodedAminoAcid.getAAs("[+42.010565]KPPK[+42.010565]YER", constants)));
		assertEquals("T,L,Y,D,E,S,C[+57],S,K[+114],E,I,Q,M[+16],A,V,L,L,K", toString(EncodedAminoAcid.getAAs("TLYDESC[+57.0214635]SK[+114.042927]EIQM[+15.994915]AVLLK", constants)));
		assertEquals("[+229],N,V,E,E,L,N,K[+229]", toString(EncodedAminoAcid.getAAs("[+229.162932]NVEELNK[+229.162932]", constants)));
		assertEquals("N,V,E,E,L,N,K", toString(EncodedAminoAcid.getAAs("NVEELNK", constants)));
	}
	
	public String toString(EncodedAminoAcid[] aas) {
		StringBuilder sb=new StringBuilder();
		for (int i = 0; i < aas.length; i++) {
			if (i>0) sb.append(",");
			sb.append(aas[i].toString());
		}
		return sb.toString();
	}
}
