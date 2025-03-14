package edu.washington.gs.maccoss.encyclopedia.filewriters.web;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2023timsTOFModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.ims.AlphaPeptDeepIMSModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt.ChronologerModel;
import junit.framework.TestCase;

public class KoinaLibraryPredictionClientTest extends TestCase {

	public void testCheckPeptide() {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		
		ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
		models.add(new Prosit2023timsTOFModel());
		models.add(new ChronologerModel());
		models.add(new AlphaPeptDeepIMSModel());
		KoinaLibraryPredictionClient client=new KoinaLibraryPredictionClient(models);

		assertFalse(checkStatus(parameters, client, "[+42.010565]VSHSELR", (byte)2));
		//123456789012345678901234567890
		//AAAAAAAAAPAAAATAPTTAATTAATAAQ
		assertTrue(checkStatus(parameters, client, "AAAAAAAAAPAAAATAPTTAATTAATAAQ", (byte)2));
		assertFalse(checkStatus(parameters, client, "[+42.010565]VNPTVFFDIAVDGEPLGR", (byte)2));

		assertFalse(checkStatus(parameters, client, "GEPLGR", (byte)2));
		assertFalse(checkStatus(parameters, client, "AAAAAAAAAPAAAATAPTTAATTAATAAQQQQ", (byte)2));
		try {
			checkStatus(parameters, client, "[+0.01]VSHSELR", (byte)2);
			fail("non standard amino acid passes");
		} catch (Exception e) {
		}
		assertFalse(checkStatus(parameters, client, "VNPT[+79.966331]VFFDIAVDGEPLGR", (byte)2));
		
		// this is intentional to support unexpected SILAC mods (there are a ton). This may not be a good call, consider for the future
		assertTrue(checkStatus(parameters, client, "VNPT[+2.010565]VFFDIAVDGEPLGR", (byte)2));
		
		try {
			checkStatus(parameters, client, "VNPT[+42.010565]VFFDIAVDGEPLGR", (byte)2);
			fail("non standard amino acid passes");
		} catch (Exception e) {
		}
	}

	public boolean checkStatus(SearchParameters parameters, KoinaLibraryPredictionClient client, String sequence,
			byte charge) {
		AminoAcidEncoding[] encoding = AminoAcidEncoding.getAAs(sequence, parameters.getAAConstants());
		boolean passes=client.checkPeptide(encoding, charge);
		return passes;
	}
}
