package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation;

import java.net.MalformedURLException;
import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.filewriters.web.PrositFragmentationPredictionModel;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class Prosit2023timsTOFModel extends PrositFragmentationPredictionModel {
	@Override
	public String getName() {
		return "Prosit 2023 timsTOF";
	}

	@Override
	public URL getURL() {
		try {
			return new URL("https://koina.wilhelmlab.org:443/v2/models/Prosit_2023_intensity_timsTOF/infer");
		} catch (MalformedURLException e) {
			throw new EncyclopediaException("Error getting Koina URL", e);
		}
	}

	@Override
	public boolean useNCE() {
		return true;
	}

}
