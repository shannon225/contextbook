package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt;

import java.net.MalformedURLException;
import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.CommonModelConstraints;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.RTPredictionModel;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class Prosit2019RTModel extends RTPredictionModel {
	@Override
	public String getName() {
		return "Prosit 2019 iRT";
	}

	@Override
	public URL getURL() {
		try {
			return new URL("https://koina.wilhelmlab.org/v2/models/Prosit_2019_irt/infer");
		} catch (MalformedURLException e) {
			throw new EncyclopediaException("Error getting Koina URL", e);
		}
	}
	
	@Override
	public boolean canModelPeptide(AminoAcidEncoding[] aas, byte precursorCharge) {
		return CommonModelConstraints.canModelPeptidePrositStandard(aas, precursorCharge);
	}

}
