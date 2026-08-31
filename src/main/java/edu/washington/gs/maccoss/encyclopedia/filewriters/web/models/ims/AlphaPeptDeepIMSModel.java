package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.ims;

import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.CommonModelConstraints;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.IMSPredictionModel;

public class AlphaPeptDeepIMSModel extends IMSPredictionModel {
	@Override
	public String getName() {
		return "AlphaPeptDeep CCS";
	}

	@Override
	public URL getURL(String baseURL) {
		return KoinaFeaturePredictionModel.inferenceURL(baseURL, "AlphaPeptDeep_ccs_generic");
	}
	
	@Override
	public boolean canModelPeptide(AminoAcidEncoding[] aas, byte precursorCharge) {
		return CommonModelConstraints.canModelPeptideAlphaPeptDeep(aas, precursorCharge);
	}

}
