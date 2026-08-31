package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt;

import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.CommonModelConstraints;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.RTPredictionModel;

public class Prosit2019RTModel extends RTPredictionModel {
	@Override
	public String getName() {
		return "Prosit 2019 iRT";
	}

	@Override
	public URL getURL(String baseURL) {
		return KoinaFeaturePredictionModel.inferenceURL(baseURL, "Prosit_2019_irt");
	}
	
	@Override
	public boolean canModelPeptide(AminoAcidEncoding[] aas, byte precursorCharge) {
		return CommonModelConstraints.canModelPeptidePrositStandard(aas, precursorCharge);
	}

}
