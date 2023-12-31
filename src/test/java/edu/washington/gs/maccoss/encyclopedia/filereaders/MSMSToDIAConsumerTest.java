package edu.washington.gs.maccoss.encyclopedia.filereaders;

import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.FloatPair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExternalResource;
import org.junit.rules.RuleChain;
import org.junit.rules.TemporaryFolder;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.io.File;
import java.util.*;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;

import static org.junit.Assert.*;

@RunWith(Parameterized.class)
public class MSMSToDIAConsumerTest {
	final int blockPrecursors;
	final int blockStripes;

	final TemporaryFolder temporaryFolder = new TemporaryFolder();

	final ConsumerRule consumerRule;

	@Rule
	public RuleChain rule;

	public MSMSToDIAConsumerTest(int blockPrecursors, int blockStripes) {
		this.blockPrecursors = blockPrecursors;
		this.blockStripes = blockStripes;

		consumerRule = new ConsumerRule();
		rule = RuleChain.outerRule(temporaryFolder)
				.around(consumerRule);
	}

	@Parameterized.Parameters
	public static Collection<Object[]> parameters() {
		return Arrays.asList(new Object[][]{
//				{100, 100}, // Outrageously slow

				// Current defaults
				{MzmlSAXToMSMSProducer.MAX_PRECURSORS_PER_BLOCK, MzmlSAXToMSMSProducer.MAX_STRIPES_PER_SCAN},

				// Original defaults
				{100, 1000}, // put this here explicitly in case the default changes

				// Larger sizes
				{100, 2000},
				{100, 5000}, // approx 2x faster than original defaults
				{1000, 5000},
				{1000, 10000},
		});
	}

	@Test
	public void testConsumer() {
		final MSMSToDIAConsumer consumer = consumerRule.getConsumer();
		assertNotNull(consumer);

		final Stopwatch stopwatch = Stopwatch.createStarted();

		final Thread thread = new Thread(consumer);
		try {
			try {
				thread.start();
				thread.join();
			} finally {
				// Stop the stopwatch immediately
				stopwatch.stop();

				// Check that everything actually exited; if not, try to clean up.
				if (thread.isAlive()) {
					thread.interrupt();
					thread.join(1000L);

					if (thread.isAlive()) {
						Logger.errorLine("Failed to join() worker thread after 1000ms; giving up!");
						fail("Thread under test never exited!");
					}
				}
			}
		} catch (InterruptedException e) {
			Logger.errorLine("Test run interrupted! Exiting without cleanup!");
			Logger.errorException(e);

			fail("Interrupted!");
		}

		assertFalse(consumer.hadError());

		// TODO: assert that .DIA contents are correct

		// Only log the time here, so we're sure we succeeded
		Logger.logLine(String.format(
				"Wrote %d blocks to .DIA in %dms (up to %d MS1 or %d MS2 per block)",
				consumerRule.numBlocks,
				stopwatch.elapsed(TimeUnit.MILLISECONDS),
				blockPrecursors,
				blockStripes
		));
	}

	private class ConsumerRule extends ExternalResource {
		private static final int NUM_SCANS = 65536;
		private static final int NUM_WINDOWS = 40;

		public static final float PRECURSOR_RANGE_LOWER = 400f;
		public static final float PRECURSOR_RANGE_UPPER = 1000f;

		private final List<FloatPair> WINDOWS = generateWindows(NUM_WINDOWS, PRECURSOR_RANGE_LOWER, PRECURSOR_RANGE_UPPER);

		private final Random random = new Random();

		private StripeFile stripeFile;

		private MSMSToDIAConsumer consumer;

		private int scanIndex = 0;

		int numBlocks;

		public final MSMSToDIAConsumer getConsumer() {
			return consumer;
		}

		@Override
		protected void before() throws Throwable {
			final File f = temporaryFolder.newFile();

			stripeFile = new StripeFile(true);
			stripeFile.openFile(f);

			Logger.logLine(String.format("Generating %d scans", NUM_SCANS));

			final Iterator<PrecursorScan> precursorStream = generatePrecursors().iterator();
			final Iterator<FragmentScan> stripeStream = generateStripes().iterator();

			final ArrayList<MSMSBlock> blocks = Lists.newArrayList();
			ArrayList<PrecursorScan> precursors = Lists.newArrayList();
			ArrayList<FragmentScan> fragments = Lists.newArrayList();
			int nPrecursors = 0;
			int nStripes = 0;
			for (int i = 0; i < NUM_SCANS; i++) {
				if (i % NUM_WINDOWS == 0) {
					precursors.add(precursorStream.next());
					nPrecursors += 1;
				}
				if (precursors.size() >= blockPrecursors) {
					blocks.add(new MSMSBlock(precursors, fragments));
					precursors = Lists.newArrayList();
					fragments = Lists.newArrayList();
				}

				fragments.add(stripeStream.next());
				nStripes += 1;

				if (fragments.size() >= blockStripes) {
					blocks.add(new MSMSBlock(precursors, fragments));
					precursors = Lists.newArrayList();
					fragments = Lists.newArrayList();
				}
			}
			if (precursors.size() + fragments.size() > 0) {
				blocks.add(new MSMSBlock(precursors, fragments));
			}

			numBlocks = blocks.size();

			Logger.logLine(String.format("Generated %d blocks of up to %d precursor or %d fragment scans", numBlocks, blockPrecursors, blockStripes));
			Logger.logLine(String.format("Average block sizes: %.02f MS1, %.02f MS2", nPrecursors / ((float) numBlocks), nStripes / ((float) numBlocks)));

			final BlockingQueue<MSMSBlock> queue = new LinkedBlockingQueue<>(numBlocks + 1);
			blocks.forEach(e -> {
				try {
					queue.put(e);
				} catch (InterruptedException ex) {
					// Won't happen -- we have the necessary capacity
					throw new IllegalStateException(ex);
				}
			});

			queue.put(MSMSBlock.POISON_BLOCK); // CRITICAL to ensure the consumer exits!

			Logger.logLine("Finished generating mock data");

			consumer = new MSMSToDIAConsumer(queue, stripeFile, SearchParameterParser.getDefaultParametersObject());
		}

		@Override
		protected void after() {
			stripeFile.close();
		}

		private Stream<PrecursorScan> generatePrecursors() {
			return Stream.generate(this::generatePrecursor);
		}

		private PrecursorScan generatePrecursor() {
			final int idx = ++scanIndex;
			final int nPeaks = 10 + random.nextInt(100);

			return new PrecursorScan(
					"spectrum" + idx,
					idx,
					random.nextFloat() * 1000f, // scan start time
					0, // fraction
					PRECURSOR_RANGE_LOWER, // isolation lower
					PRECURSOR_RANGE_UPPER, // isolation upper
					null, // ion inject time
					generateMasses(nPeaks, 400f, 1000f),
					generateIntensities(nPeaks),
					Optional.empty()
			);
		}

		private Stream<FragmentScan> generateStripes() {
			return Stream.generate(() -> this.generateStripe(WINDOWS));
		}

		private FragmentScan generateStripe(List<FloatPair> windows) {
			final int idx = ++scanIndex;
			final int nPeaks = 10 + random.nextInt(100);

			final FloatPair window = windows.get(idx % windows.size());

			return new FragmentScan(
					"spectrum" + idx,
					null,
					idx,
					random.nextFloat() * 1000f, // scan start time
					0, // fraction
					null, // ion inject time
					window.getOne(),
					window.getTwo(),
					generateMasses(nPeaks, 400f, 1000f),
					generateIntensities(nPeaks),
					Optional.empty()
			);
		}

		private double[] generateMasses(int n, double min, double max) {
			return DoubleStream.generate(() -> min + random.nextDouble() * (max - min))
					.limit(n)
					.toArray();
		}

		private float[] generateIntensities(int n) {
			return General.toFloatArray(
					DoubleStream.generate(() -> random.nextDouble() * 1e9)
							.limit(n)
							.toArray()
			);
		}

		private List<FloatPair> generateWindows(int numWindows, float min, float max) {
			final List<FloatPair> result = Lists.newArrayListWithCapacity(numWindows);

			final float width = (max - min) / ((float) numWindows);

			for (int i = 0; i < numWindows; i++) {
				result.add(new FloatPair(min + i * width, min + (i + 1) * width));
			}

			return result;
		}
	}
}
