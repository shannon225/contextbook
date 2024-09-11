package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation;

import java.net.MalformedURLException;
import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.filewriters.web.PrositFragmentationPredictionModel;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class Prosit2020CIDModel extends PrositFragmentationPredictionModel {
	@Override
	public String getName() {
		return "Prosit 2020 CID";
	}

	@Override
	public URL getURL(String baseURL) {
		try {
			return new URL(baseURL+"v2/models/Prosit_2020_intensity_CID/infer");
		} catch (MalformedURLException e) {
			throw new EncyclopediaException("Error getting Koina URL", e);
		}
	}

	@Override
	public boolean useNCE() {
		return false;
	}

}
