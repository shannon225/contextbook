package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import junit.framework.TestCase;

public class EncodedAminoAcidTest extends TestCase {
	public void testEncoding() {
		AminoAcidConstants constants=new AminoAcidConstants();
		assertEquals(",Q[-17],A,H,L,C[+57],V,L,A,S,N,C[+57],D,E,P,M,Y,V,K,", toString(AminoAcidEncoding.getAAs("Q[-17.026549]AHLC[+57.021464]VLASNC[+57.021464]DEPMYVK", constants)));
		assertEquals("[+42],S[+80],G,S,S,S,V,A,A,M[+16],K,", toString(AminoAcidEncoding.getAAs("[+42.010565]S[+79.966331]GSSSVAAM[+15.994915]K", constants)));
		assertEquals("[+42],S,L,L,D,G,L,A,S,S[+80],P,R,A,P,L,Q,S,S,K,", toString(AminoAcidEncoding.getAAs("[+42.010565]SLLDGLASS[+79.966331]PRAPLQSSK", constants)));
		assertEquals(",N,D,I,K[+42],L,A,A,K[+42],L,I,H,T,L,D,D,R,", toString(AminoAcidEncoding.getAAs("NDIK[+42.010565]LAAK[+42.010565]LIHTLDDR", constants)));
		assertEquals(",K[+42],P,P,K[+42],Y,E,R,", toString(AminoAcidEncoding.getAAs("K[+42.010565]PPK[+42.010565]YER", constants)));
		assertEquals("[+42],K,P,P,K[+42],Y,E,R,", toString(AminoAcidEncoding.getAAs("[+42.010565]KPPK[+42.010565]YER", constants)));
		assertEquals(",T,L,Y,D,E,S,C[+57],S,K[+114],E,I,Q,M[+16],A,V,L,L,K,", toString(AminoAcidEncoding.getAAs("TLYDESC[+57.0214635]SK[+114.042927]EIQM[+15.994915]AVLLK", constants)));
		assertEquals("[+229],N,V,E,E,L,N,K[+229],", toString(AminoAcidEncoding.getAAs("[+229.162932]NVEELNK[+229.162932]", constants)));
		assertEquals(",N,V,E,E,L,N,K,", toString(AminoAcidEncoding.getAAs("NVEELNK", constants)));
	}
	public void testPeptideModSeq() {
		AminoAcidConstants constants=new AminoAcidConstants();
		assertEquals("Q[-17.026549]AHLC[+57.0214635]VLASNC[+57.0214635]DEPMYVK", AminoAcidEncoding.getPeptideModSeq(AminoAcidEncoding.getAAs("Q[-17.026549]AHLC[+57.021464]VLASNC[+57.021464]DEPMYVK", constants), constants));
		assertEquals("[+42.010565]S[+79.966331]GSSSVAAM[+15.994915]K", AminoAcidEncoding.getPeptideModSeq(AminoAcidEncoding.getAAs("[+42.010565]S[+79.966331]GSSSVAAM[+15.994915]K", constants), constants));
		assertEquals("[+42.010565]SLLDGLASS[+79.966331]PRAPLQSSK", AminoAcidEncoding.getPeptideModSeq(AminoAcidEncoding.getAAs("[+42.010565]SLLDGLASS[+79.966331]PRAPLQSSK", constants), constants));
		assertEquals("NDIK[+42.010565]LAAK[+42.010565]LIHTLDDR", AminoAcidEncoding.getPeptideModSeq(AminoAcidEncoding.getAAs("NDIK[+42.010565]LAAK[+42.010565]LIHTLDDR", constants), constants));
		assertEquals("K[+42.010565]PPK[+42.010565]YER", AminoAcidEncoding.getPeptideModSeq(AminoAcidEncoding.getAAs("K[+42.010565]PPK[+42.010565]YER", constants), constants));
		assertEquals("[+42.010565]KPPK[+42.010565]YER", AminoAcidEncoding.getPeptideModSeq(AminoAcidEncoding.getAAs("[+42.010565]KPPK[+42.010565]YER", constants), constants));
		assertEquals("TLYDESC[+57.0214635]SK[+114.042927]EIQM[+15.994915]AVLLK", AminoAcidEncoding.getPeptideModSeq(AminoAcidEncoding.getAAs("TLYDESC[+57.0214635]SK[+114.042927]EIQM[+15.994915]AVLLK", constants), constants));
		//assertEquals("N[+229.162932]VEELNK[+229.162932]", AminoAcidEncoding.toString(AminoAcidEncoding.getAAs("[+229.162932]NVEELNK[+229.162932]", constants), constants));
		assertEquals("NVEELNK", AminoAcidEncoding.getPeptideModSeq(AminoAcidEncoding.getAAs("NVEELNK", constants), constants));
	}
	
	public String toString(AminoAcidEncoding[] aas) {
		StringBuilder sb=new StringBuilder();
		for (int i = 0; i < aas.length; i++) {
			if (i>0) sb.append(",");
			sb.append(aas[i].toString());
		}
		return sb.toString();
	}
}
