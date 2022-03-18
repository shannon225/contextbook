package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExternalResource;
import org.junit.rules.RuleChain;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.util.concurrent.LinkedBlockingQueue;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;

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

		final Thread thread = new Thread(consumer);
		try {
			try {
				thread.start();
				thread.join();
			} finally {
				if (thread.isAlive()) {
					thread.interrupt();
				}
				thread.join();
			}
		} catch (InterruptedException e) {
			Logger.errorLine("Test run interrupted! Exiting without cleanup!");
			Logger.errorException(e);
		}

		assertFalse(consumer.hadError());
	}

	private static class ConsumerRule extends ExternalResource {
		private final TemporaryFolder temporaryFolder;

		private MSMSToDIAConsumer consumer;

		public ConsumerRule(TemporaryFolder temporaryFolder) {
			this.temporaryFolder = temporaryFolder;
		}

		@Override
		protected void before() throws Throwable {
			final File f = temporaryFolder.newFile();

			final StripeFile stripeFile = new StripeFile(true);
			stripeFile.openFile(f);

			final LinkedBlockingQueue<MSMSBlock> queue = new LinkedBlockingQueue<>();

			// TODO: Populate queue with chunks prior to starting the test

			queue.put(MSMSBlock.POISON_BLOCK); // CRITICAL to ensure the consumer exits!

			consumer = new MSMSToDIAConsumer(queue, stripeFile, SearchParameterParser.getDefaultParametersObject());
		}

		public final MSMSToDIAConsumer getConsumer() {
			return consumer;
		}
	}
}
