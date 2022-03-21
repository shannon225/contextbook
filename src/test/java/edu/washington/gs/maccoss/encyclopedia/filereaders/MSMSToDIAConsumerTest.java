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

import java.io.File;
import java.util.List;
import java.util.Random;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;

import static org.junit.Assert.*;

public class MSMSToDIAConsumerTest {
	TemporaryFolder temporaryFolder = new TemporaryFolder();

	ConsumerRule consumerRule = new ConsumerRule(temporaryFolder);

	@Rule
	public RuleChain rule = RuleChain.outerRule(temporaryFolder)
			.around(consumerRule);

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
				"Wrote %d blocks to .DIA in %dms using NUMBER_OF_STRIPES_AT_ONCE = %d",
				ConsumerRule.NUM_BLOCKS,
				stopwatch.elapsed(TimeUnit.MILLISECONDS),
				StripeFile.NUMBER_OF_STRIPES_AT_ONCE
		));
	}

	private static class ConsumerRule extends ExternalResource {
		private static final int NUM_BLOCKS = 65536 / (MzmlSAXToMSMSProducer.MAX_PRECURSORS_PER_BLOCK + MzmlSAXToMSMSProducer.MAX_STRIPES_PER_SCAN);
		private static final int NUM_WINDOWS = 40;

		public static final float PRECURSOR_RANGE_LOWER = 400f;
		public static final float PRECURSOR_RANGE_UPPER = 1000f;

		private static final List<FloatPair> WINDOWS = generateWindows(NUM_WINDOWS, PRECURSOR_RANGE_LOWER, PRECURSOR_RANGE_UPPER);

		private final Random random = new Random();

		private final TemporaryFolder temporaryFolder;

		private MSMSToDIAConsumer consumer;

		private int scanIndex = 0;

		public ConsumerRule(TemporaryFolder temporaryFolder) {
			this.temporaryFolder = temporaryFolder;
		}

		public final MSMSToDIAConsumer getConsumer() {
			return consumer;
		}

		@Override
		protected void before() throws Throwable {
			final File f = temporaryFolder.newFile();

			final StripeFile stripeFile = new StripeFile(true);
			stripeFile.openFile(f);

			final BlockingQueue<MSMSBlock> queue = new LinkedBlockingQueue<>(NUM_BLOCKS + 1);

			Logger.logLine(String.format("Generating %d blocks of %d precursor and %d fragment scans", NUM_BLOCKS, MzmlSAXToMSMSProducer.MAX_PRECURSORS_PER_BLOCK, MzmlSAXToMSMSProducer.MAX_STRIPES_PER_SCAN));

			generateMsMsBlocks().limit(NUM_BLOCKS).forEach(e -> {
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

		private Stream<MSMSBlock> generateMsMsBlocks() {
			return Stream.generate(this::generateMsMsBlock);
		}

		private MSMSBlock generateMsMsBlock() {
			List<PrecursorScan> precursors = generatePrecursors()
					.limit(MzmlSAXToMSMSProducer.MAX_PRECURSORS_PER_BLOCK)
					.collect(Collectors.toList());

			List<FragmentScan> stripes = generateStripes()
					.limit(MzmlSAXToMSMSProducer.MAX_STRIPES_PER_SCAN)
					.collect(Collectors.toList());

			return new MSMSBlock(
					precursors,
					stripes
			);
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
					generateIntensities(nPeaks)
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
					generateIntensities(nPeaks)
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

		private static List<FloatPair> generateWindows(int numWindows, float min, float max) {
			final List<FloatPair> result = Lists.newArrayListWithCapacity(numWindows);

			final float width = (max - min) / ((float) numWindows);

			for (int i = 0; i < numWindows; i++) {
				result.add(new FloatPair(min + i * width, min + (i + 1) * width));
			}

			return result;
		}
	}
}
