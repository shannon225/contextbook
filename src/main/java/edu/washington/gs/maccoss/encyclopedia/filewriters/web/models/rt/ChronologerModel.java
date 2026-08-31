package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt;

import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.CommonModelConstraints;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.RTPredictionModel;

public class ChronologerModel extends RTPredictionModel {
	@Override
	public String getName() {
		return "Chronologer %ACN";
	}
	
	@Override
	public String getDataTypeName() {
		return "rt";
	}

	@Override
	public URL getURL(String baseURL) {
		return KoinaFeaturePredictionModel.inferenceURL(baseURL, "Chronologer_RT");
	}
	
	@Override
	public boolean canModelPeptide(AminoAcidEncoding[] aas, byte precursorCharge) {
		return CommonModelConstraints.canModelChronologer(aas, precursorCharge);
	}

}
