package edu.washington.gs.maccoss.encyclopedia.utils.io;

import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;

import static org.junit.Assert.assertNotNull;

public class TableParserTest {
	private static final long TIMEOUT = 5000L;

	private Path tmp;

	@Before
	public void setUp() throws Exception {
		tmp = Files.createTempFile("test_", ".tsv");
	}

	@After
	public void tearDown() throws Exception {
		FileUtils.deleteQuietly(tmp.toFile());
		tmp = null;
	}

	@Test(timeout = TIMEOUT)
	public void testSimpleMuscle() throws Exception {
		TableParser.parseTSV(
				tmp.toFile(),
				simpleMuscle()
		);
	}

	@Test(timeout = TIMEOUT)
	public void testNoSuchFile() throws Exception {

	}

	@Test(timeout = TIMEOUT)
	public void testEmptyFile() throws Exception {

	}

	@Test(timeout = TIMEOUT)
	public void testErrorProducing() throws Exception {

	}

	@Test(timeout = TIMEOUT)
	public void testErrorConsuming() throws Exception {

	}

	static TableParserMuscle simpleMuscle() {
		return new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				assertNotNull(row);
			}

			@Override
			public void cleanup() { }
		};
	}
}