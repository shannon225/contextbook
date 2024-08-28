package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation;

import java.net.MalformedURLException;
import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.filewriters.web.PrositFragmentationPredictionModel;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class Prosit2020HCDModel extends PrositFragmentationPredictionModel {
	@Override
	public String getName() {
		return "Prosit 2020 HCD";
	}

	@Override
	public URL getURL() {
		try {
			return new URL("https://koina.wilhelmlab.org/v2/models/Prosit_2020_intensity_HCD/infer");
		} catch (MalformedURLException e) {
			throw new EncyclopediaException("Error getting Koina URL", e);
		}
	}

	@Override
	public boolean useNCE() {
		return true;
	}

}
