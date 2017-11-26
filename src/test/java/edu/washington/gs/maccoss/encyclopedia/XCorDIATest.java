package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;
import java.util.zip.DataFormatException;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.ScoringResultsToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import gnu.trove.map.hash.TCharDoubleHashMap;
import junit.framework.AssertionFailedError;
import junit.framework.TestCase;

/**
 * @author Seth Just
 * @since 2:20 PM 10/16/17
 */
public class XCorDIATest extends TestCase {
	/**
	 * This test attempts to run XCorDIA to see if we improperly compute peptide
	 * or fragment masses at any point. The issue is that with synthetic data it's
	 * hard to make Percolator NOT fail, so we work around that by overriding
	 * methods to process results as they're produced (I tried setting up a
	 * TeeConsumer but it was too difficult to work with the threading to ensure
	 * it was joined at the right time).
	 *
	 * The test fails if no cysteine-containing peptides are encountered by our
	 * overridden method, if any of the assertions in that method fails, or if
	 * Percolator does not encounter the expected error. Note that the test case
	 * fails or succeeds based on flags, while the "true" cause of the failure is
	 * written to the console when it's encountered.
	 *
	 * If the sanity-check exception in {@link PeptideUtils#getMasses(String, AminoAcidConstants)}
	 * is disabled, this test (currently) fails because of strange fragments in entries
	 * that don't exist in the correct.
	 *
	 * If the sanity-check exception is enabled, then no cysteine-containing peptides
	 * make it to our overridden method, and the error is written to the logs.
	 */
	public void testFixedModCalculations() throws Exception {
		PecanSearchParameters parameters=PecanParameterParser.getDefaultParametersObject();
		System.out.println("STARTING PARAMS: "+parameters.getAAConstants().getFixedModString());
		FastaEntry entry=new FastaEntry("APEPTIDEKACPEPTIDECKMARECYSPEPTIDESK");
		ArrayList<FastaPeptideEntry> seqs=parameters.getEnzyme().digestProtein(entry, parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages(), parameters.getAAConstants());

		final File fastaFile = File.createTempFile("test_", ".fasta");
		fastaFile.deleteOnExit();

		final File diaFile = File.createTempFile("test_", ".dia");
		diaFile.delete(); // ensure it doesn't exist

		final PercolatorExecutionData percolatorFiles = XCorDIAJobData.getPercolatorExecutionData(diaFile, fastaFile, parameters);

		percolatorFiles.getInputTSV().deleteOnExit();
		percolatorFiles.getProteinDecoyFile().deleteOnExit();
		percolatorFiles.getProteinOutputFile().deleteOnExit();
		percolatorFiles.getPeptideOutputFile().deleteOnExit();
		percolatorFiles.getPeptideDecoyFile().deleteOnExit();

		final AtomicBoolean processedResult = new AtomicBoolean(false);
		final AtomicBoolean failed = new AtomicBoolean(false);

		XCorDIAJobData jobData = new XCorDIAJobData(
				Optional.of(seqs),
				diaFile,
				new FakeStripeFile(diaFile),
				fastaFile,
				percolatorFiles,
				new XCorDIAOneScoringFactory(new PecanSearchParameters(
						new AminoAcidConstants(), // includes C+57 ONLY
						FragmentationType.YONLY,
						new MassTolerance(10),
						new MassTolerance(10),
						DigestionEnzyme.getEnzyme("trypsin"),
						false
				)) {
					@Override
					public PeptideScoringResultsConsumer getResultsConsumer(File outputFile, BlockingQueue<PeptideScoringResult> resultsQueue, StripeFileInterface diaFile) {
						final ScoringResultsToTSVConsumer resultsConsumer = (ScoringResultsToTSVConsumer) super.getResultsConsumer(outputFile, resultsQueue, diaFile);
						return new ScoringResultsToTSVConsumer(outputFile, diaFile, resultsConsumer.getScoreNames(), resultsQueue, 1) {

							// Contains NO fixed mods
							private final AminoAcidConstants NO_MODS = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());

							@Override
							protected void processResult(PeptideScoringResult result) {
								try {
									final LibraryEntry entry = result.getEntry();
									if (entry.getPeptideSeq().contains("C")) {
										processedResult.set(true);
										
										assertTrue(
												"Cysteine-containing peptide had no +57 in modSeq: " + entry.getPeptideModSeq(),
												entry.getPeptideModSeq().contains("C[+57")
										);

										assertEquals(
												"Cysteine-containing peptide had incorrect m/z!",
												NO_MODS.getChargedMass(entry.getPeptideModSeq(), entry.getPrecursorCharge()),
												entry.getPrecursorMZ(),
												0.0001d
										);

										// Get fragments based on the fixed mods in the modSeq
										final FragmentIon[] fragmentIons = new FragmentationModel(entry.getPeptideModSeq(), NO_MODS)
												.getPrimaryIonObjects(FragmentationType.YONLY, (byte)1); // XCorr model only has +1 ions
										final double[] fragMasses = Arrays.stream(fragmentIons)
												.mapToDouble(ion -> ion.mass)
												.toArray();

										double[] modelMasses=entry.getMassArray();
										Arrays.stream(fragMasses)
												.forEach(frag -> { // for each mass in the entry's mass array...
													assertTrue( // assert that a matching mass is in the model's fragments
															"Fragment ion " + frag + " was not found in fragmentation model for "
																	+ entry.getPeptideModSeq() + "+" + entry.getPrecursorCharge()
																	+ " (model=" + Arrays.toString(fragMasses) + ")"
															,
															Arrays.stream(modelMasses)
																	.anyMatch(modelMass ->
																			// check within a wide tolerance
																			new MassTolerance(100)
																					.compareTo(frag, modelMass) == 0
																	)||frag>SparseXCorrCalculator.biggestFragmentMass
													);
												});
									}
								} catch (AssertionFailedError err) {
									err.printStackTrace();
									failed.set(true);
								}

								super.processResult(result);
							}
						};
					}
				}
		);

		try {
			XCorDIA.runPie(new EmptyProgressIndicator(), jobData);
		} catch (EncyclopediaException e) {
			assertTrue(
					"Unexpected failure: \"" + e.getMessage() + "\"",
					e.getMessage().contains("Percolator exited with non-zero status: 139")
			);
		}

		assertTrue("Did not process any results for cysteine-containing peptides!", processedResult.get());
		assertFalse("Encountered failure during processing!", failed.get());
	}

	private FastaPeptideEntry toEntry(String pepSeq) {
		return new FastaPeptideEntry(pepSeq);
	}

	private static class FakeStripeFile implements StripeFileInterface {
		private final File file;

		private FakeStripeFile(File file) {
			this.file = file;
		}

		@Override
		public Map<Range, Float> getRanges() {
			return ImmutableMap.of(new Range(100f, 1200f), 10f);
		}

		@Override
		public void openFile(File userFile) throws IOException, SQLException {

		}

		@Override
		public ArrayList<PrecursorScan> getPrecursors(float minRT, float maxRT) throws IOException, SQLException, DataFormatException {
			return Lists.newArrayList();
		}

		@Override
		public ArrayList<Stripe> getStripes(double targetMz, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException {
			return getStripes(new Range((float) targetMz - 1f, (float) targetMz + 1f), minRT, maxRT, sqrt);
		}

		@Override
		public ArrayList<Stripe> getStripes(Range targetMzRange, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException {
			double[] masses = new double[] {
					 147.11285d, // y1 of ACPEPTIDECK (no fixed mod)
					 329.12785d, // b3 of ACPEPTIDECK (includes fixed mod)
					 531.27082d // b4 of MORE... (no fixed mod)
			};

			final float[] intens = new float[masses.length];
			Arrays.fill(intens, 1f);

			return Lists.newArrayList(
					new Stripe("0", "0", 0, 0f, 400f, 700f, masses, intens),
					new Stripe("1", "1", 1, 1f, 400f, 700f, masses, intens),
					new Stripe("2", "2", 2, 2f, 400f, 700f, masses, intens),
					new Stripe("3", "3", 3, 3f, 400f, 700f, masses, intens),
					new Stripe("4", "4", 4, 4f, 400f, 700f, masses, intens),
					new Stripe("5", "5", 5, 5f, 400f, 700f, masses, intens),
					new Stripe("6", "6", 6, 6f, 400f, 700f, masses, intens),
					new Stripe("7", "7", 7, 7f, 400f, 700f, masses, intens)
			);
		}

		@Override
		public float getTIC() throws IOException, SQLException {
			return 1f;
		}

		@Override
		public float getGradientLength() throws IOException, SQLException {
			return 100f;
		}

		@Override
		public void close() {
		}

		@Override
		public boolean isOpen() {
			return true;
		}

		@Override
		public File getFile() {
			return file;
		}

		@Override
		public String getOriginalFileName() {
			return file.getName();
		}
	}
}