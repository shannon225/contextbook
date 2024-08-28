package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation;

import java.net.MalformedURLException;
import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.filewriters.web.IMSPredictionModel;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class IM2DeepModel extends IMSPredictionModel {
	@Override
	public String getName() {
		return "IM2Deep CCS";
	}

	@Override
	public URL getURL() {
		try {
			return new URL("https://koina.wilhelmlab.org:443/v2/models/IM2Deep/infer");
		} catch (MalformedURLException e) {
			throw new EncyclopediaException("Error getting Koina URL", e);
		}
	}

}
