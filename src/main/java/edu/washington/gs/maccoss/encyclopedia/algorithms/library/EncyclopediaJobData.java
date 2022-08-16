package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
import edu.washington.gs.maccoss.encyclopedia.filereaders.*;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Map;
import java.util.zip.DataFormatException;

public class EncyclopediaJobData extends QuantitativeSearchJobData implements LibrarySearchJobData {
	public static final String LOG_FILE_SUFFIX=".log";
	public static final String DECOY_PROTEIN_FILE_SUFFIX=".encyclopedia.protein_decoy.txt";
	public static final String OUTPUT_PROTEIN_FILE_SUFFIX=".encyclopedia.protein.txt";
	public static final String DECOY_FILE_SUFFIX=".encyclopedia.decoy.txt";
	public static final String OUTPUT_FILE_SUFFIX=".encyclopedia.txt";
	public static final String FEATURE_FILE_SUFFIX=".features.txt";

	private final LibraryInterface library;
	private final LibraryScoringFactory taskFactory;

	public EncyclopediaJobData(File diaFile, File fastaFile, LibraryInterface library, LibraryScoringFactory taskFactory) {
		this(diaFile, null, getPercolatorExecutionData(diaFile, fastaFile, taskFactory.getParameters()), taskFactory.getParameters(), ProgramType.getGlobalVersion().toString(), library, taskFactory);
	}

	public EncyclopediaJobData(File diaFile, File fastaFile, LibraryInterface library, File outputFile, LibraryScoringFactory taskFactory) {
		this(diaFile, null, getPercolatorExecutionData(outputFile, fastaFile, taskFactory.getParameters()), taskFactory.getParameters(), ProgramType.getGlobalVersion().toString(), library, taskFactory);
	}

	public EncyclopediaJobData(File diaFile, PercolatorExecutionData percolatorFiles, SearchParameters parameters, String version, LibraryInterface library, LibraryScoringFactory taskFactory) {
		this(diaFile, null, percolatorFiles, parameters, version, library, taskFactory);
	}

	public EncyclopediaJobData(File diaFile, StripeFileInterface diaFileReader, PercolatorExecutionData percolatorFiles, SearchParameters parameters, String version, LibraryInterface library, LibraryScoringFactory taskFactory) {
		super(diaFile, diaFileReader, percolatorFiles, parameters, version);

		this.library = library;
		this.taskFactory = taskFactory;
	}

	public static PercolatorExecutionData getPercolatorExecutionData(File referenceFileLocation, File fastaFile, SearchParameters parameters) {
		return new PercolatorExecutionData(new File(getPrefixFromOutput(referenceFileLocation) + FEATURE_FILE_SUFFIX), fastaFile,
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_FILE_SUFFIX), 
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_PROTEIN_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_PROTEIN_FILE_SUFFIX), parameters);
	}

	static String getPrefixFromOutput(File outputFile) {
		final String absolutePath = outputFile.getAbsolutePath();

		if (absolutePath.endsWith(OUTPUT_FILE_SUFFIX)) {
			return absolutePath.substring(0, absolutePath.length() - OUTPUT_FILE_SUFFIX.length());
		} else {
			return absolutePath;
		}
	}

	public EncyclopediaJobData updateTaskFactory(LibraryScoringFactory taskFactory) {
		return new EncyclopediaJobData(getDiaFile(), diaFileReader, getPercolatorFiles(), getParameters(), getVersion(), getLibrary(), taskFactory);
	}

	public LibraryInterface getLibrary() {
		return library;
	}

	public LibraryScoringFactory getTaskFactory() {
		return taskFactory;
	}

	public File getResultLibrary() {
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile());
		return new File(absolutePath + LibraryFile.ELIB);
	}

	@Override
	public String getSearchType() {
		return "EncyclopeDIA";
	}
	
	@Override
	public String getPrimaryScoreName() {
		return taskFactory.getPrimaryScoreName();
	}

	/**
	 * Return an {@code EncyclopeDIAJobData} instance that for {@code diaFile} that will
	 * function even if the .DIA file doesn't exist. Only useful when the file's results
	 * are present but the .DIA file doesn't. This allows SearchToBLIB's {@code -alignOnly}
	 * option to run without access to raw data, which is useful e.g. for large experiments
	 * where collecting all the raw files is costly.
	 * <p>
	 * The returned instance will return {@code true} from {@link #hasBeenRun()} even if
	 * the .DIA doesn't exist. It will also return a {@code StripeFileInterface} that provides
	 * the appropriate original file name (read from the file's results ELIB; see
	 * {@link StripeFileInterface#getOriginalFileName()}) but otherwise throws on attempts to
	 * read the .DIA's contents.
	 */
	public static EncyclopediaJobData getDummyFor(File diaFile, File fastaFile, LibraryInterface library, LibraryScoringFactory factory) {
		return new DummyEncyclopediaJobData(diaFile, fastaFile, library, factory);
	}

	/**
	 * Special class overriding key methods to allow processing results for jobs
	 * even if the DIA file doesn't exist. This allows SearchToBLIB's {@code -alignOnly}
	 * option to run without access to raw data, which is useful e.g. for large
	 * experiments where collecting all the raw files is costly.
	 * <p>
	 * Instances will return {@code true} from {@link #hasBeenRun()} even if
	 * ths .DIA doesn't exist. It will also return a {@code StripeFileInterface} that provides
	 * the appropriate original file name (read from the file's results ELIB; see
	 * {@link StripeFileInterface#getOriginalFileName()}) but otherwise throws on attempts to
	 * read the .DIA's contents.
	 */
	private static class DummyEncyclopediaJobData extends EncyclopediaJobData {
		private String originalFileName = null;

		private DummyEncyclopediaJobData(File diaFile, File fastaFile, LibraryInterface library, LibraryScoringFactory factory) {
			super(diaFile, fastaFile, library, factory);

			if (Files.exists(diaFile.toPath())) {
				Logger.errorLine("Creating a dummy job datum for " + diaFile.getName() + " even though it exists!");
			}
		}

		@Override
		public boolean hasBeenRun() {
			final PercolatorExecutionData percolatorFiles = getPercolatorFiles();
			if (!percolatorFiles.getInputTSV().exists()) {
				Logger.errorLine("Missing feature file: " + percolatorFiles.getInputTSV().getName());
				return false;
			}
			if (!percolatorFiles.getPeptideOutputFile().exists()) {
				Logger.errorLine("Missing output file: " + percolatorFiles.getPeptideOutputFile().getName());
				return false;
			}
			if (!getResultLibrary().exists()) {
				Logger.errorLine("Missing output library: " + getResultLibrary().getName());
				return false;
			}
			return true;
		}

		@Override
		public StripeFileInterface getDiaFileReader() {
			return new StripeFileInterface() {
				@Override
				public Map<Range, WindowData> getRanges() {
					throw new UnsupportedOperationException("File not found: " + getDiaFile().getAbsolutePath());
				}

				@Override
				public void openFile(File userFile) throws IOException, SQLException {
					throw new UnsupportedOperationException();
				}

				@Override
				public ArrayList<PrecursorScan> getPrecursors(float minRT, float maxRT) throws IOException, SQLException, DataFormatException {
					throw new UnsupportedOperationException("File not found: " + getDiaFile().getAbsolutePath());
				}

				@Override
				public ArrayList<FragmentScan> getStripes(double targetMz, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException {
					throw new UnsupportedOperationException("File not found: " + getDiaFile().getAbsolutePath());
				}

				@Override
				public ArrayList<FragmentScan> getStripes(Range targetMzRange, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException {
					throw new UnsupportedOperationException("File not found: " + getDiaFile().getAbsolutePath());
				}

				@Override
				public float getTIC() throws IOException, SQLException {
					throw new UnsupportedOperationException("File not found: " + getDiaFile().getAbsolutePath());
				}

				@Override
				public float getGradientLength() throws IOException, SQLException {
					throw new UnsupportedOperationException("File not found: " + getDiaFile().getAbsolutePath());
				}

				@Override
				public void close() {
					// no-op
				}

				@Override
				public boolean isOpen() {
					return false;
				}

				@Override
				public File getFile() {
					return getDiaFile();
				}

				@Override
				public String getOriginalFileName() {
					if (null != originalFileName) {
						return originalFileName;
					}

					synchronized (this) {
						if (null == originalFileName) {
							// Workaround: if the DIA file is missing we can't read the
							// original file name, which likely has an .mzML extension.
							// Instead, get the name used in this job's results ELIB.
							try (Connection c = new SQLFile() {}.getConnection(getResultLibrary())) {
								try (Statement s = c.createStatement()) {
									try (ResultSet rs = s.executeQuery(
											"SELECT sourcefile" +
													" FROM entries" +
													" LIMIT 1;"
									)) {
										if (rs.next()) {
											return rs.getString(1);
										} else {
											throw new SQLException("No entries in results ELIB!");
										}
									}
								}
							} catch (IOException | SQLException e) {
								Logger.errorLine("Unable to read from results ELIB for job " + getDiaFile().getName());
								Logger.errorException(e);
							}

							originalFileName = getDiaFile().getName();
						}
					}

					return originalFileName;
				}
			};
		}
	}
}