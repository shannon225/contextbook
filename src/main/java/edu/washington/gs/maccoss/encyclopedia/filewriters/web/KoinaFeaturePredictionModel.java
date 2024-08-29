package edu.washington.gs.maccoss.encyclopedia.filewriters.web;

import java.net.URL;
import java.util.List;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;

public interface KoinaFeaturePredictionModel {
	public String getName();
	public URL getURL();
	public void updatePeptides(List<KoinaPrecursor> peptides);
	public boolean canModelPeptide(AminoAcidEncoding[] aas, byte precursorCharge);
}
