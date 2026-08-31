package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation;

import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.PrositFragmentationPredictionModel;

public class Prosit2023timsTOFModel extends PrositFragmentationPredictionModel {
	@Override
	public String getName() {
		return "Prosit 2023 timsTOF";
	}

	@Override
	public URL getURL(String baseURL) {
		return KoinaFeaturePredictionModel.inferenceURL(baseURL, "Prosit_2023_intensity_timsTOF");
	}

	@Override
	public boolean useNCE() {
		return true;
	}

}
