package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt;

import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.CommonModelConstraints;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.RTPredictionModel;

public class Prosit2020TMTRTModel extends RTPredictionModel {
	@Override
	public String getName() {
		return "Prosit 2020 TMT iRT";
	}

	@Override
	public URL getURL(String baseURL) {
		return KoinaFeaturePredictionModel.inferenceURL(baseURL, "Prosit_2020_irt_TMT");
	}
	
	@Override
	public boolean canModelPeptide(AminoAcidEncoding[] aas, byte precursorCharge) {
		return CommonModelConstraints.canModelPeptidePrositTMT(aas, precursorCharge);
	}

}
