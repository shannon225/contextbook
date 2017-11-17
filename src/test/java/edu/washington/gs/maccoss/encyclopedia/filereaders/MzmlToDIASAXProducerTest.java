package edu.washington.gs.maccoss.encyclopedia.filereaders;

import com.google.common.collect.ImmutableMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
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
		// TODO Decode
	}

}