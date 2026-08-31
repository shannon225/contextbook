package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt;

import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.CommonModelConstraints;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.RTPredictionModel;

public class DeepLCHelaRTModel extends RTPredictionModel {
	@Override
	public String getName() {
		return "DeepLC HeLa iRT";
	}

	@Override
	public URL getURL(String baseURL) {
		return KoinaFeaturePredictionModel.inferenceURL(baseURL, "Deeplc_hela_hf");
	}
	
	@Override
	public boolean canModelPeptide(AminoAcidEncoding[] aas, byte precursorCharge) {
		return CommonModelConstraints.canModelPeptideDeepLC(aas, precursorCharge);
	}

}
