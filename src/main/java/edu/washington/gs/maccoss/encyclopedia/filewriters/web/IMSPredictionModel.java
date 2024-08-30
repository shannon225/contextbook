package edu.washington.gs.maccoss.encyclopedia.filewriters.web;

import java.util.List;

public abstract class IMSPredictionModel extends SingleValuePredictionModel {

	@Override
	public String getDataTypeName() {
		return "ccs";
	}
	
	@Override
	public boolean useCharge() {
		return true;
	}
	
	@Override
	public boolean useNCE() {
		return false;
	}
	
	@Override
	public void updateResults(List<KoinaPrecursor> peptides, float[] values) {
		assert(peptides.size()==values.length);
		for (int i = 0; i < peptides.size(); i++) {
			peptides.get(i).setIMS(values[i]);
		}
	}
}
