package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation;

import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.PrositFragmentationPredictionModel;

public class Prosit2020CIDModel extends PrositFragmentationPredictionModel {
	@Override
	public String getName() {
		return "Prosit 2020 CID";
	}

	@Override
	public URL getURL(String baseURL) {
		return KoinaFeaturePredictionModel.inferenceURL(baseURL, "Prosit_2020_intensity_CID");
	}

	@Override
	public boolean useNCE() {
		return false;
	}

}
