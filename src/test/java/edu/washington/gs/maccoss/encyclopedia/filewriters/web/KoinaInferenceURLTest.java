package edu.washington.gs.maccoss.encyclopedia.filewriters.web;

import java.net.URL;

import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2020HCDModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.ims.IM2DeepIMSModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt.Prosit2019RTModel;
import junit.framework.TestCase;

public class KoinaInferenceURLTest extends TestCase {
	private static final String EXPECTED="https://koina.wilhelmlab.org:443/v2/models/Prosit_2019_irt/infer";

	public void testTrailingSlashIsOptional() {
		assertEquals(EXPECTED,
				KoinaFeaturePredictionModel.inferenceURL("https://koina.wilhelmlab.org:443/", "Prosit_2019_irt").toString());
		assertEquals(EXPECTED,
				KoinaFeaturePredictionModel.inferenceURL("https://koina.wilhelmlab.org:443", "Prosit_2019_irt").toString());
	}

	public void testEveryModelTypeBuildsTheSameWayWithoutASlash() {
		String base="https://koina.wilhelmlab.org:443";

		URL fragmentation=new Prosit2020HCDModel().getURL(base);
		URL ims=new IM2DeepIMSModel().getURL(base);
		URL rt=new Prosit2019RTModel().getURL(base);

		assertEquals(base+"/v2/models/Prosit_2020_intensity_HCD/infer", fragmentation.toString());
		assertEquals(base+"/v2/models/IM2Deep/infer", ims.toString());
		assertEquals(base+"/v2/models/Prosit_2019_irt/infer", rt.toString());
	}

	public void testEveryRegisteredModelProducesAUsableURL() {
		for (KoinaFeaturePredictionModel model : KoinaFeaturePredictionModel.getFragmentationModels()) {
			assertURLIsWellFormed(model);
		}
		for (KoinaFeaturePredictionModel model : KoinaFeaturePredictionModel.getIMSModels()) {
			assertURLIsWellFormed(model);
		}
		for (KoinaFeaturePredictionModel model : KoinaFeaturePredictionModel.getRTModels()) {
			assertURLIsWellFormed(model);
		}
	}

	private static void assertURLIsWellFormed(KoinaFeaturePredictionModel model) {
		URL url=model.getURL("https://koina.wilhelmlab.org:443");
		String message=model.getName()+" built ["+url+"]";
		assertTrue(message, url.getPath().startsWith("/v2/models/"));
		assertTrue(message, url.getPath().endsWith("/infer"));
		assertEquals(message, "koina.wilhelmlab.org", url.getHost());
		assertEquals(message, 443, url.getPort());
	}
}
