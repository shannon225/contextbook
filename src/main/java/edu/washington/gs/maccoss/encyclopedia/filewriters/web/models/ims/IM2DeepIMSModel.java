package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.ims;

import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.CommonModelConstraints;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.IMSPredictionModel;

public class IM2DeepIMSModel extends IMSPredictionModel {
	@Override
	public String getName() {
		return "IM2Deep CCS";
	}

	@Override
	public URL getURL(String baseURL) {
		return KoinaFeaturePredictionModel.inferenceURL(baseURL, "IM2Deep");
	}
	
	@Override
	public boolean canModelPeptide(AminoAcidEncoding[] aas, byte precursorCharge) {
		return CommonModelConstraints.canModelPeptideDeepLC(aas, precursorCharge);
	}
}
