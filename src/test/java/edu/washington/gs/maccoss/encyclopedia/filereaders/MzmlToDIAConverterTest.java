package edu.washington.gs.maccoss.encyclopedia.filereaders;

import com.google.common.base.Stopwatch;
import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.concurrent.TimeUnit;

import static org.junit.Assert.*;

public class MzmlToDIAConverterTest {
	private static final int[] QUEUE_CAPACITIES = { 1, 4, 16, 64 };

	@Rule
	public TemporaryFolder temporaryFolder = new TemporaryFolder();

	@Test
	public void benchMzmlToDIAConverter() throws Throwable {
		final ImmutableList.Builder<Path> inputPathsBuilder = ImmutableList.builder();

		try (InputStream is = getClass().getResourceAsStream("mzml_to_dia_bench_inputs.txt")) {
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

		final boolean origPrint = Logger.PRINT_TO_SCREEN;
		try {
			Logger.PRINT_TO_SCREEN = false;

			System.out.println("file\tsize (MiB)\tqueueCapacity\ttime (ms)");
			inputPaths.forEach(this::benchConversion);
		} finally {
			Logger.PRINT_TO_SCREEN = origPrint;
		}
	}

	private void benchConversion(Path mzmlPath) {
		for (int cap : QUEUE_CAPACITIES) {
			benchConversion(mzmlPath, cap);
		}
	}

	private void benchConversion(Path mzmlPath, int queueCapacity) {
		final Stopwatch stopwatch = Stopwatch.createStarted();
		try {
			MzmlToDIAConverter.convertSAX(
					mzmlPath.toFile(),
					temporaryFolder.newFile(),
					SearchParameterParser.getDefaultParametersObject(),
					false, // crashes with true!!
					queueCapacity
			);
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		} finally {
			stopwatch.stop();
		}

		try {
			System.out.println(String.format(
					"%s\t%.02f\t%d\t%d",
					mzmlPath.getFileName().toString(),
					Files.size(mzmlPath) / ((float) 1024*1024),
					queueCapacity,
					stopwatch.elapsed(TimeUnit.MILLISECONDS))
			);
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}
	}
}