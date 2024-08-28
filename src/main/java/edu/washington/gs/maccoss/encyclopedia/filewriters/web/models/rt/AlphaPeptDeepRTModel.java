package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt;

import java.net.MalformedURLException;
import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.filewriters.web.RTPredictionModel;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class AlphaPeptDeepRTModel extends RTPredictionModel {
	@Override
	public String getName() {
		return "AlphaPeptDeep iRT";
	}

	@Override
	public URL getURL() {
		try {
			return new URL("https://koina.wilhelmlab.org:443/v2/models/AlphaPeptDeep_rt_generic/infer");
		} catch (MalformedURLException e) {
			throw new EncyclopediaException("Error getting Koina URL", e);
		}
	}

}
