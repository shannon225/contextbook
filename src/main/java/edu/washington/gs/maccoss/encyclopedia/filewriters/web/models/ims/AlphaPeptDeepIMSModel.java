package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.ims;

import java.net.MalformedURLException;
import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.CommonModelConstraints;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.IMSPredictionModel;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class AlphaPeptDeepIMSModel extends IMSPredictionModel {
	@Override
	public String getName() {
		return "AlphaPeptDeep CCS";
	}

	@Override
	public URL getURL() {
		try {
			return new URL("https://koina.wilhelmlab.org:443/v2/models/AlphaPeptDeep_ccs_generic/infer");
		} catch (MalformedURLException e) {
			throw new EncyclopediaException("Error getting Koina URL", e);
		}
	}
	
	@Override
	public boolean canModelPeptide(AminoAcidEncoding[] aas, byte precursorCharge) {
		return CommonModelConstraints.canModelPeptideAlphaPeptDeep(aas, precursorCharge);
	}

}
