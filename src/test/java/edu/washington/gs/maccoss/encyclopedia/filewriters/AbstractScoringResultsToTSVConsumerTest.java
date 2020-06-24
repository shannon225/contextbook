package edu.washington.gs.maccoss.encyclopedia.filewriters;

import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import junit.framework.TestCase;

import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.nio.file.Path;
import java.nio.file.Paths;

public class AbstractScoringResultsToTSVConsumerTest extends TestCase {
	public void testSortingPerformanceAndTuning() throws Exception {
		final ImmutableList.Builder<Path> inputPathsBuilder = ImmutableList.builder();

		try (InputStream is = getClass().getResourceAsStream("tsv_sorting_inputs.txt")) {
			try (LineNumberReader reader = new LineNumberReader(new InputStreamReader(is))) {
				String line;
				while (null != (line = reader.readLine())) {
					inputPathsBuilder.add(Paths.get(line));
				}
			}
		}

		final ImmutableList<Path> inputPaths = inputPathsBuilder.build();

		Logger.logLine(String.format("Read %d paths from test resource file.", inputPaths.size()));

		inputPaths.forEach(this::testSortTuning);
	}

	private void testSortTuning(Path path) {
		fail("Not yet implemented"); //TODO
	}
}