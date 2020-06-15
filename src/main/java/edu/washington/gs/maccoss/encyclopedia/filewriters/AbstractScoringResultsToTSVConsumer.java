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

		Logger.logLine("Sorting results into " + this.outputFile.getAbsolutePath());

		try {
			final Serializer<CSVRecord> serializer = Serializer.csv(
					CSVFormat
							.DEFAULT
							.withDelimiter('\t')
							.withRecordSeparator(getLineSeparator())
							.withFirstRecordAsHeader(),
					StandardCharsets.UTF_8
			);

			// compare rows on column zero (ID) using lexicographic string ordering
			final Comparator<CSVRecord> comparator = Comparator.comparing(
					record -> record.get(0),
					Comparator.naturalOrder()
			);

			Sorter
					.serializer(serializer)
					.comparator(comparator)
					.input(this.tmpFile)
					.output(this.outputFile)
					.tempDirectory(this.outputFile.getParentFile()) // always sort in the target directory
					.maxItemsPerFile(100000) // 100k lines per file; this controls memory usage
					.sort();
		} catch (UncheckedIOException exception) {
			Logger.errorLine("Caught IO exception sorting TSV output; failing!");
			Logger.errorException(exception);
			throw exception;
		} finally {
			// unconditionally remove the unsorted temp file
			System.out.println("Removing temp file " + this.tmpFile.getAbsolutePath());
			FileUtils.deleteQuietly(this.tmpFile);
		}
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
