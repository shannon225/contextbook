package edu.washington.gs.maccoss.encyclopedia.utils.io;

import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import static org.junit.Assert.*;

public class TableParserTest {
	private static final long TIMEOUT = 5000L;

	private Path tmp;

	private AtomicInteger rowCount;
	private AtomicBoolean didCleanup;
	private TableParserMuscle muscle = new TableParserMuscle() {
		@Override
		public void processRow(Map<String, String> row) {
			assertNotNull(row);
			rowCount.incrementAndGet();
		}

		@Override
		public void cleanup() {
			assertTrue(didCleanup.compareAndSet(false, true));
		}
	};

	@Before
	public void setUp() throws Exception {
		tmp = Files.createTempFile("test_", ".tsv");
		rowCount = new AtomicInteger(0);
		didCleanup = new AtomicBoolean(false);
	}

	@After
	public void tearDown() throws Exception {
		FileUtils.deleteQuietly(tmp.toFile());
		tmp = null;
		rowCount = null;
		didCleanup = null;
	}

	@Test(timeout = TIMEOUT)
	public void testSimpleMuscle() throws Exception {
		//TODO: set up file contents

		TableParser.parseTSV(
				tmp.toFile(),
				muscle
		);

//		assertEquals(, rowCount.get());
		assertTrue(didCleanup.get());
	}

	@Test(timeout = TIMEOUT, expected = FileNotFoundException.class)
	public void testNoSuchFile() throws Throwable {
		Files.delete(tmp); // throws if not deleted

		assertFalse(Files.exists(tmp));

		try {
			TableParser.parseTSV(
					tmp.toFile(),
					muscle
			);
		} catch (EncyclopediaException e) {
			// We'll get any IO exception wrapped in an EncyclopediaException,
			// unwrapping it allows us to expect a more specific exception
			// from this test case.
			throw e.getCause() == null ? e : e.getCause();
		}

		assertEquals(0, rowCount.get());
		assertFalse(didCleanup.get());
	}

	@Test(timeout = TIMEOUT)
	public void testEmptyFile() throws Exception {
		Files.write(tmp, new byte[0]);

		assertTrue(Files.exists(tmp));

		TableParser.parseTSV(
				tmp.toFile(),
				muscle
		);

		assertEquals(0, rowCount.get());
		assertFalse(didCleanup.get());
	}

	@Test(timeout = TIMEOUT, expected = Exception.class) //TODO: assert more specific exception type
	public void testErrorProducing() throws Exception {
		Files.write(
				tmp,
				ImmutableList.of(
						"A\tB\tC",
						"a\tb\tc",
						"a\tb\tc"
				),
				StandardCharsets.UTF_8
		);

		assertTrue(Files.exists(tmp));

		TableParser.parseTable(
				muscle,
				new TableParserProducer(new LinkedBlockingQueue<>(), tmp.toFile(), "\t", 1) {
					@Override
					public void run() {
						throw new CustomException();
					}
				}
		);

		// We can process no lines or one line before the error,
		// but it shouldn't be possible to process more.
		final int rowCount = this.rowCount.get();
		assertEquals("Processed unexpected number of rows", 0, rowCount);

		assertFalse(didCleanup.get());
	}

	@Test(timeout = TIMEOUT, expected = Exception.class) //TODO: assert more specific exception type
	public void testErrorConsuming() throws Exception {
		Files.write(
				tmp,
				ImmutableList.of(
						"A\tB\tC",
						"a\tb\tc"
				),
				StandardCharsets.UTF_8
		);

		assertTrue(Files.exists(tmp));

		TableParser.parseTSV(
				tmp.toFile(),
				new TableParserMuscle() {
					@Override
					public void processRow(Map<String, String> row) {
						throw new CustomException();
					}

					@Override
					public void cleanup() {
						assertTrue(didCleanup.compareAndSet(false, true));
					}
				}
		);

		assertTrue(didCleanup.get());
	}

	static class CustomException extends RuntimeException { }
}