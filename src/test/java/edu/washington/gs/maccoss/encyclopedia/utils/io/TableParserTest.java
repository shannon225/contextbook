package edu.washington.gs.maccoss.encyclopedia.utils.io;

import junit.framework.TestCase;
import org.apache.commons.io.FileUtils;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;

public class TableParserTest extends TestCase {
	private Path tmp;

	public void setUp() throws Exception {
		tmp = Files.createTempFile("test_", ".tsv");
	}

	public void tearDown() throws Exception {
		FileUtils.deleteQuietly(tmp.toFile());
		tmp = null;
	}

	public void testSimpleMuscle() throws Exception {
		TableParser.parseTSV(
				tmp.toFile(),
				simpleMuscle()
		);
	}

	public void testNoSuchFile() throws Exception {

	}

	public void testEmptyFile() throws Exception {

	}

	public void testErrorProducing() throws Exception {

	}

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