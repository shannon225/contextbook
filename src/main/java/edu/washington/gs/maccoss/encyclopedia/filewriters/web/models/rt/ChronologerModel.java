package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt;

import java.net.MalformedURLException;
import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.CommonModelConstraints;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.RTPredictionModel;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class ChronologerModel extends RTPredictionModel {
	@Override
	public String getName() {
		return "Chronologer %ACN";
	}
	
	@Override
	public String getDataTypeName() {
		return "rt";
	}

	@Override
	public URL getURL(String baseURL) {
		try {
			return new URL(baseURL+"v2/models/Chronologer_RT/infer");
		} catch (MalformedURLException e) {
			throw new EncyclopediaException("Error getting Koina URL", e);
		}
	}
	
	@Override
	public boolean canModelPeptide(AminoAcidEncoding[] aas, byte precursorCharge) {
		return CommonModelConstraints.canModelChronologer(aas, precursorCharge);
	}

}
