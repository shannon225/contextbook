package edu.washington.gs.maccoss.encyclopedia.utils.threading;

import java.util.Random;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;
import org.hamcrest.core.IsEqual;
import org.hamcrest.core.IsNot;
import org.junit.Assume;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class ExternalExecutorTest {
	public static void main(String[] args) throws Exception {
		ExternalExecutor e=new ExternalExecutor(new String[] {"cat", "/Users/searleb/Documents/data/dbs/UP000005640_9606.fasta"});
		BlockingQueue<OutputMessage> result=e.start();

		int count=0;
		while (!e.isFinished()||!result.isEmpty()) {
			OutputMessage data=result.take();
			if (data.getMessage().startsWith(">")) {
				count++;
				if (count%1000==0) System.out.println(count+", "+result.size()+", "+e.isFinished());
			}
		}
		System.out.println(count+", "+result.size());
		System.out.println("FINISHED!");
	}

	@Test
	public void testExternalOutputRaceCondition() throws Exception {
		// In code prior to 2.0.0-SNAPSHOT there's a potential race condition if the
		// output queue is ever drained entirely after the process has exited, even
		// if there's still unprocessed output from the process! This can conceivably
		// happen to any process, but realistically will only affect very short-lived
		// processes, or output that's written very close to when the process exits.
		//
		// The root cause of this bug is that ExternalExecutor#isFinished() returns
		// true if the process has exited but the output-handling threads are still
		// running!

		// Skip the test on Windows
		Assume.assumeThat(OSDetector.getOS(), new IsNot<>(new IsEqual<>(OSDetector.OS.WINDOWS)));

		// Build a random test string
		final int n = new Random().nextInt(100);
		final StringBuilder sb = new StringBuilder();
		for (int i = n; i --> 0;) {
			sb.append(i);
			if (i > 0) {
				sb.append("\n");
			}
		}

		final ExternalExecutor e = new ExternalExecutor(new String[] {"echo", sb.toString()});

		final BlockingQueue<?> q = e.start();

		int count = 0;

		// If we uncomment this line this test case almost always succeeds, as there's enough
		// time for at least one message to be enqueued.
//		Thread.sleep(1);

		// Match existing usages of ExternalExecutor -- loop until finished and no output
		while (!e.isFinished() || !q.isEmpty()) {
			if (!q.isEmpty()) {
				q.take();
				count++;
			}
		}

		assertEquals(n, count);
	}
}
