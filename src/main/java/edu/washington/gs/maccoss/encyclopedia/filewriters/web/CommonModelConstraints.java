package edu.washington.gs.maccoss.encyclopedia.filewriters.web;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;

public class CommonModelConstraints {

	public static boolean canModelPeptidePrositStandard(AminoAcidEncoding[] aas, byte precursorCharge) {
		if (precursorCharge<1||precursorCharge>6) {
			return false;
		}
		if (aas.length<7||aas.length>30) return false;
		
		for (int i = 0; i < aas.length; i++) {
			if (aas[i].getIndex()>22) return false;
		}
		return true;
	}

	public static boolean canModelPeptidePrositTMT(AminoAcidEncoding[] aas, byte precursorCharge) {
		if (precursorCharge<1||precursorCharge>6) {
			return false;
		}
		if (aas.length>30) return false;
		
		// Prosit-TMT requires n-term TMT, but not K-TMT!
		if (aas[0]!=AminoAcidEncoding.nTMT10) return false;
		
		for (int i = 1; i < aas.length; i++) {
			if (aas[i].getIndex()>22) {
				if (aas[i]!=AminoAcidEncoding.KTMT10) { 
					return false;
				}
			}
		}
		return true;
	}
	public static boolean canModelPeptideAlphaPeptDeep(AminoAcidEncoding[] aas, byte precursorCharge) {
		if (precursorCharge<2||precursorCharge>4) {
			return false;
		}
		if (aas.length<7||aas.length>35) return false;
		
		for (int i = 0; i < aas.length; i++) {
			if (aas[i].getIndex()>22) {
				if (aas[i]!=AminoAcidEncoding.nAc) { 
					return false;
				}
			}
		}
		return true;
	}
	public static boolean canModelPeptideDeepLC(AminoAcidEncoding[] aas, byte precursorCharge) {
		return true; // TODO is this really the case?!
	}
}
