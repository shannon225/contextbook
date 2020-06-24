package edu.washington.gs.maccoss.encyclopedia.filewriters;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import junit.framework.TestCase;

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;

public class AbstractScoringResultsToTSVConsumerTest extends TestCase {
	private static final ImmutableList<Integer> SORT_SIZES = ImmutableList.of(30000, 100000, 300000, 1000000, 3000000);

	public void testSortingPerformanceAndTuning() throws Exception {
		final ImmutableList.Builder<Path> inputPathsBuilder = ImmutableList.builder();

		try (InputStream is = getClass().getResourceAsStream("tsv_sorting_inputs.txt")) {
			try (LineNumberReader reader = new LineNumberReader(new InputStreamReader(is))) {
				String line;
				while (null != (line = reader.readLine())) {
					// ignore if blank or first non-whitespace is '#'
					// https://regexper.com/#%5Cs*%28%23.*%29%3F
					if (!line.matches("\\s*(#.*)?")) {
						inputPathsBuilder.add(Paths.get(line));
					}
				}
			}
		}

		final ImmutableList<Path> inputPaths = inputPathsBuilder.build();

		Logger.logLine(String.format("Read %d paths from test resource file.", inputPaths.size()));
		Logger.timelessLogLine("");

		Logger.timelessLogLine("file\tsize\t"+ Joiner.on(" (ms)\t").join(SORT_SIZES) + " (ms)");

		inputPaths.forEach(this::testSortTuning);
	}

	private void testSortTuning(Path path) {
		final File file = path.toFile();

		assertTrue(file.exists());
		assertTrue(file.canRead());

		Logger.log(path.getFileName()+"\t"+file.length());

		final File output;
		try {
			output = File.createTempFile(
					"output_", ".tsv",
					file.getParentFile() // in same folder (so on same volume usually)
			);
			output.deleteOnExit();
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}

		SORT_SIZES.forEach(i -> {
			final long nanos = AbstractScoringResultsToTSVConsumer.doFileSort(i, file, output, System.lineSeparator());
			Logger.log("\t" + nanos / 1e6);

			output.delete();
		});

		Logger.timelessLogLine("");
	}
}