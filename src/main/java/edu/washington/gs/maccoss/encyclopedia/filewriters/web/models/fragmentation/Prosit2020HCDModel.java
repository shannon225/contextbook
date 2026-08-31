package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation;

import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.PrositFragmentationPredictionModel;

public class Prosit2020HCDModel extends PrositFragmentationPredictionModel {
	@Override
	public String getName() {
		return "Prosit 2020 HCD";
	}

	@Override
	public URL getURL(String baseURL) {
		return KoinaFeaturePredictionModel.inferenceURL(baseURL, "Prosit_2020_intensity_HCD");
	}

	@Override
	public boolean useNCE() {
		return true;
	}

}
