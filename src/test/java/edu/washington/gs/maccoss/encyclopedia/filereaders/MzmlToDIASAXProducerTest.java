package edu.washington.gs.maccoss.encyclopedia.filereaders;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultimap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.mzml.InstrumentComponent;
import edu.washington.gs.maccoss.encyclopedia.filereaders.mzml.InstrumentId;
import edu.washington.gs.maccoss.encyclopedia.filereaders.mzml.InstrumentMapTranscoder;
import junit.framework.TestCase;
import org.mockito.Mockito;

import java.io.File;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class MzmlToDIASAXProducerTest extends TestCase {

	interface FileAndExpectedSoftwareVersions {
		File getTestMzMLFile() throws URISyntaxException;
		Map<String, String> getExpectedSoftwareAccessionIdToVersion();
	}

	private FileAndExpectedSoftwareVersions getTestResource() {
		return new FileAndExpectedSoftwareVersions() {
			@Override
			public File getTestMzMLFile() throws URISyntaxException {
				URL resource = getClass().getClassLoader().getResource("truncated-mzml.mzml");
				if (resource == null) {
					throw new NullPointerException("Could not find resource!");
				} else {
					return new File(resource.toURI());
				}
			}

			@Override
			public Map<String, String> getExpectedSoftwareAccessionIdToVersion() {
				return ImmutableMap.of(
						"MS:1000615", "3.0.9987",
						"MS:1000551", "unknown");
			}
		};
	}

	public void testGetSoftwareAccessionIdToVersion() throws Exception {

		FileAndExpectedSoftwareVersions harness = getTestResource();

		File mzMLFile = harness.getTestMzMLFile();
		Map<String, String> expected = harness.getExpectedSoftwareAccessionIdToVersion();

		File diaFileSaveDestination = new File(mzMLFile.getParentFile(), "test-dia" + StripeFile.DIA_EXTENSION);
		diaFileSaveDestination.deleteOnExit();

		MzmlToDIAConverter.convertSAX(mzMLFile, diaFileSaveDestination, Mockito.mock(SearchParameters.class), false);
		StripeFile f = new StripeFile();
		f.openFile(diaFileSaveDestination);
		HashMap<String, String> metadata = f.getMetadata();

		Set<String> softwareVersions = metadata.entrySet()
				.stream()
				.filter(entry -> entry.getKey().startsWith(StripeFile.SOFTWARE_VERSION_PREFIX))
				.map(entry -> entry.getKey().substring(StripeFile.SOFTWARE_VERSION_PREFIX.length()) + "=" + entry.getValue())
				.collect(Collectors.toSet());

		assertTrue(softwareVersions.size() == expected.size());
		expected.forEach((key, value) -> {
			assertTrue(softwareVersions.contains(key + "=" + value));
		});

	}

	public void testGetInstrumentConfiguration() throws Exception {
		FileAndExpectedSoftwareVersions harness = getTestResource();
		File mzMLFile = harness.getTestMzMLFile();
		File diaFileSaveDestination = new File(mzMLFile.getParentFile(), "test-dia" + StripeFile.DIA_EXTENSION);
		MzmlToDIAConverter.convertSAX(mzMLFile, diaFileSaveDestination, Mockito.mock(SearchParameters.class), false);
		StripeFile f = new StripeFile();
		f.openFile(diaFileSaveDestination);
		HashMap<String, String> metadata = f.getMetadata();
		String encodedInstrumentConfigurations = metadata.get(StripeFile.INSTRUMENT_CONFIGURATIONS);
		ImmutableMultimap<InstrumentId, InstrumentComponent> decodedInstrumentMap = InstrumentMapTranscoder.decode(encodedInstrumentConfigurations);
		assertTrue(decodedInstrumentMap.asMap().size() == 1);
		decodedInstrumentMap.asMap().forEach((entry, value) -> {
			assertEquals("IC1", entry.instrumentConfigurationId);
			assertEquals("MS:1002533", entry.accession);
			assertEquals("TripleTOF 6600", entry.name);
			value.forEach(component -> {
				switch (component.order) {
					case 1:
						assertEquals("MS", component.cvRef);
						assertEquals("MS:1000073", component.accessionId);
						assertEquals("electrospray ionization", component.name);
						break;
					case 2:
						assertEquals("MS", component.cvRef);
						assertEquals("MS:1000081", component.accessionId);
						assertEquals("quadrupole", component.name);
						break;
					case 3:
						assertEquals("MS", component.cvRef);
						assertEquals("MS:1000081", component.accessionId);
						assertEquals("quadrupole", component.name);
						break;
					case 4:
						assertEquals("MS", component.cvRef);
						assertEquals("MS:1000084", component.accessionId);
						assertEquals("time-of-flight", component.name);
						break;
					case 5:
						assertEquals("MS", component.cvRef);
						assertEquals("MS:1000253", component.accessionId);
						assertEquals("electron multiplier", component.name);
						break;
					default:
						fail();
				}
			});
		});
	}

}