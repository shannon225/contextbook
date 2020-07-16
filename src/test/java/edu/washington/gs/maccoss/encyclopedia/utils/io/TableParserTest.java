package edu.washington.gs.maccoss.encyclopedia.utils.io;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.Map;

import static org.junit.Assert.*;

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
		//TODO: set up file contents

		TableParser.parseTSV(
				tmp.toFile(),
				simpleMuscle()
		);
	}

	@Test(timeout = TIMEOUT, expected = FileNotFoundException.class)
	public void testNoSuchFile() throws Throwable {
		Files.delete(tmp); // throws if not deleted

		assertFalse(Files.exists(tmp));

		try {
			TableParser.parseTSV(
					tmp.toFile(),
					simpleMuscle()
			);
		} catch (EncyclopediaException e) {
			// We'll get any IO exception wrapped in an EncyclopediaException,
			// unwrapping it allows us to expect a more specific exception
			// from this test case.
			throw e.getCause() == null ? e : e.getCause();
		}
	}

	@Test(timeout = TIMEOUT)
	public void testEmptyFile() throws Exception {
		Files.write(tmp, new byte[0], StandardOpenOption.TRUNCATE_EXISTING);

		assertTrue(Files.exists(tmp));

		TableParser.parseTSV(
				tmp.toFile(),
				simpleMuscle()
		);
	}

	@Test(timeout = TIMEOUT, expected = Exception.class) //TODO: assert more specific exception type
	public void testErrorProducing() throws Exception {

	}

	@Test(timeout = TIMEOUT, expected = Exception.class) //TODO: assert more specific exception type
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