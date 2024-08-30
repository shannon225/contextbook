package edu.washington.gs.maccoss.encyclopedia.filewriters.web;

import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;

public class KoinaPrecursor {
	private final AminoAcidEncoding[] aas;
	private final float nce;
	private final byte charge;
	private volatile float iRT=Float.NEGATIVE_INFINITY;
	private volatile float IMS=Float.NEGATIVE_INFINITY;
	private volatile float[] intensities=null;
	private volatile double[] mzs=null;
	HashSet<String> accessions=new HashSet<String>();
	
	public KoinaPrecursor(AminoAcidEncoding[] aas, float nce, byte charge) {
		this.aas=aas;
		this.nce = nce;
		this.charge = charge;
	}
	
	public String getKoinaSequence() {
		return AminoAcidEncoding.getUnimodSeq(aas);
	}
	public HashSet<String> getAccessions() {
		return accessions;
	}
	public byte getCharge() {
		return charge;
	}
	public float getNCE() {
		return nce;
	}
	
	public void addAccession(String accession) {
		accessions.add(accession);
	}

	public void addAccessions(HashSet<String> multipleaccessions) {
		accessions.addAll(multipleaccessions);
	}
	
	public void setIMS(float IMS) {
		this.IMS = IMS;
	}
	
	public void setiRT(float iRT) {
		this.iRT = iRT;
	}
	
	public void setIntensities(float[] intensities) {
		this.intensities = intensities;
	}
	
	public void setMzs(double[] mzs) {
		this.mzs = mzs;
	}
	
	@Override
	public int hashCode() {
		return getKoinaSequence().hashCode()+16807*Float.hashCode(nce)+1771561*Byte.hashCode(charge);
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj==null||!(obj instanceof KoinaPrecursor)) {
			return false;
		}
		KoinaPrecursor k=(KoinaPrecursor)obj;
		if (aas.length!=k.aas.length) return false;
		for (int i = 0; i < aas.length; i++) {
			if (aas[i]!=k.aas[i]) return false;
		}
		if (nce!=k.nce) {
			return false;
		}
		if (charge!=k.charge) {
			return false;
		}
		return true;
	}
	
	public AnnotatedLibraryEntry toEntry(AminoAcidConstants constants, SearchParameters params) {
		String peptideModSeq=AminoAcidEncoding.getPeptideModSeq(aas, constants);
		double precursorMZ=constants.getChargedMass(peptideModSeq, charge);
		
		Optional<Float> maybeIMS;
		if (Float.isInfinite(IMS)) {
			maybeIMS=Optional.empty();
		} else {
			maybeIMS=Optional.of(IMS);
		}
		
		float rtInSec;
		if (Float.isInfinite(iRT)) {
			rtInSec=0.0f;
		} else {
			rtInSec=iRT*60f;
		}
		if (mzs==null) {
			mzs=new double[0];
			intensities=new float[0];
		}
		LibraryEntry e=new LibraryEntry("Koina", accessions, precursorMZ, charge, peptideModSeq, 1, rtInSec,
				0.0f, mzs, intensities, maybeIMS, constants);
		
		return new AnnotatedLibraryEntry(e, params);
	}
}
