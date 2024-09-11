package edu.washington.gs.maccoss.encyclopedia.filewriters.web;

import java.util.List;

public abstract class RTPredictionModel extends SingleValuePredictionModel {

	@Override
	public String getDataTypeName() {
		return "irt";
	}
	
	@Override
	public boolean useCharge() {
		return false;
	}
	
	@Override
	public boolean useNCE() {
		return false;
	}
	
	@Override
	public String getModelType() {
		return KoinaFeaturePredictionModel.RT_TYPE;
	}
	
	@Override
	public void updateResults(List<KoinaPrecursor> peptides, float[] values) {
		assert(peptides.size()==values.length);
		for (int i = 0; i < peptides.size(); i++) {
			peptides.get(i).setiRT(values[i]);
		}
	}
}
