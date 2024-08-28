package edu.washington.gs.maccoss.encyclopedia.filewriters.web;

import java.net.URL;
import java.util.List;

public interface KoinaFeaturePredictionModel {
	public String getName();
	public URL getURL();
	public void updatePeptides(List<KoinaPrecursor> peptides);
}
