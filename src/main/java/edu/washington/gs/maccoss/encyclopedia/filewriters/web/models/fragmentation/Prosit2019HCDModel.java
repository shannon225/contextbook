package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation;

import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.PrositFragmentationPredictionModel;

public class Prosit2019HCDModel extends PrositFragmentationPredictionModel {
	@Override
	public String getName() {
		return "Prosit 2019 HCD";
	}

	@Override
	public URL getURL(String baseURL) {
		return KoinaFeaturePredictionModel.inferenceURL(baseURL, "Prosit_2019_intensity");
	}

	@Override
	public boolean useNCE() {
		return true;
	}

}
