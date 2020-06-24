package edu.washington.gs.maccoss.encyclopedia.filewriters;

import com.github.davidmoten.bigsorter.Serializer;
import com.github.davidmoten.bigsorter.Sorter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.io.FileUtils;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.Comparator;
import java.util.concurrent.BlockingQueue;

/**
 * Base class for a consumer that writes all results (with a header) to a temp file,
 * then sorts those files (using string comparison of the first column) into the destination
 * file when {@code close()} is called.
 */
public abstract class AbstractScoringResultsToTSVConsumer implements PeptideScoringResultsConsumer {
	protected final OSDetector.OS os = OSDetector.getOS();

	protected final File outputFile;
	protected final File tmpFile;
	protected final StripeFileInterface diaFile;
	protected final BlockingQueue<PeptideScoringResult> resultsQueue;
	protected final PrintWriter writer;

	protected volatile int numberProcessed = 0;

	public AbstractScoringResultsToTSVConsumer(File outputFile, StripeFileInterface diaFile, BlockingQueue<PeptideScoringResult> resultsQueue) {
		this.outputFile = outputFile;
		this.tmpFile = new File(outputFile.getAbsolutePath() + ".unsorted");
		this.diaFile = diaFile;
		this.resultsQueue = resultsQueue;

		try {
			writer = new PrintWriter(this.tmpFile, "UTF-8");
			System.out.println("Constructing writer for " + this.tmpFile.getAbsolutePath());
		} catch (FileNotFoundException e) {
			throw new EncyclopediaException("Error setting up output file: " + this.tmpFile.getAbsolutePath(), e);
		} catch (UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error setting up output file: " + this.tmpFile.getAbsolutePath(), e);
		}
	}

	@Override
	public final BlockingQueue<PeptideScoringResult> getResultsQueue() {
		return resultsQueue;
	}

	@Override
	public final void close() {
		writer.flush();
		writer.close();

		try {
			doFileSort(tmpFile, outputFile, getLineSeparator());
		} catch (UncheckedIOException exception) {
			Logger.errorLine("Caught IO exception sorting TSV output; failing!");
			Logger.errorException(exception);
			throw exception;
		} finally {
			// unconditionally remove the unsorted temp file
			System.out.println("Removing temp file " + tmpFile.getAbsolutePath());
			FileUtils.deleteQuietly(tmpFile);
		}
	}

	static void doFileSort(File inputFile, File outputFile, String lineSeparator) throws UncheckedIOException {
		doFileSort(100000, inputFile, outputFile, lineSeparator); // 100k lines per file; this controls memory usage
	}

	static void doFileSort(int maxItems, File inputFile, File outputFile, String lineSeparator) throws UncheckedIOException {
		Logger.logLine("Sorting results into " + outputFile.getAbsolutePath());

		final Serializer<CSVRecord> serializer = Serializer.csv(
				CSVFormat
						.DEFAULT
						.withDelimiter('\t')
						.withRecordSeparator(lineSeparator)
						.withFirstRecordAsHeader(),
				StandardCharsets.UTF_8
		);

		// compare rows on column zero (ID) using lexicographic string ordering
		final Comparator<CSVRecord> comparator = Comparator.comparing(
				record -> record.get(0),
				Comparator.naturalOrder()
		);

		final long start = System.nanoTime();

		Sorter
				.serializer(serializer)
				.comparator(comparator)
				.input(inputFile)
				.output(outputFile)
				.tempDirectory(outputFile.getParentFile()) // always sort in the target directory
				.maxItemsPerFile(maxItems) // this controls memory usage
				.sort();

		final long end = System.nanoTime();
		Logger.logLine(String.format("Sorted feature file in %.02f seconds wall clock time", (end - start) / 1e9));
	}

	@Override
	public final int getNumberProcessed() {
		return numberProcessed;
	}

	@Override
	public abstract void run();

	private String getLineSeparator() {
		switch (os) {
			case MAC:
				return "\n";
			default:
				return System.lineSeparator();
		}
	}
}
