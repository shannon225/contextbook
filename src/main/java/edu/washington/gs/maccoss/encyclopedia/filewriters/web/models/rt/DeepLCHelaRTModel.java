package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt;

import java.net.MalformedURLException;
import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.filewriters.web.RTPredictionModel;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class DeepLCHelaRTModel extends RTPredictionModel {
	@Override
	public String getName() {
		return "DeepLC HeLa iRT";
	}

	@Override
	public URL getURL() {
		try {
			return new URL("https://koina.wilhelmlab.org:443/v2/models/Deeplc_hela_hf/infer");
		} catch (MalformedURLException e) {
			throw new EncyclopediaException("Error getting Koina URL", e);
		}
	}

}
