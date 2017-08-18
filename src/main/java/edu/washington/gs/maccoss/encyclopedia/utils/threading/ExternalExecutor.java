package edu.washington.gs.maccoss.encyclopedia.utils.threading;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;

public class ExternalExecutor {
	private final String[] cmdAndArgs;
	private volatile Process p;

	public ExternalExecutor(String[] cmdAndArgs) {
		this.cmdAndArgs=cmdAndArgs;
	}

	public BlockingQueue<OutputMessage> start() throws IOException {
		final BlockingQueue<OutputMessage> queue=new LinkedBlockingQueue<OutputMessage>();

		ImmutableList.Builder<String> cmdBuilder = ImmutableList.builder();

		final OS os=OSDetector.getOS();
		if (OS.WINDOWS==os) {
			cmdBuilder.add("cmd").add("/c");
		}

		cmdBuilder.addAll(Arrays.asList(cmdAndArgs));

		final ImmutableList<String> cmd = cmdBuilder.build();

		Logger.logLine("Executing ["+ Joiner.on(" ").join(cmd)+"]");

		p = new ProcessBuilder(cmd).start();

		// Note: both these threads should terminate by exiting the loop once the process completes
		new Thread(new Runnable() {
			public void run() {
				BufferedReader input=new BufferedReader(new InputStreamReader(p.getInputStream()));
				String line=null;

				try {
					while ((line=input.readLine())!=null) {
						queue.add(new OutputMessage(line, true));
					}
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}).start();

		new Thread(new Runnable() {
			public void run() {
				BufferedReader input=new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line=null;

				try {
					while ((line=input.readLine())!=null) {
						queue.add(new OutputMessage(line, false));
					}
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}).start();

		return queue;
	}

	public boolean isFinished() {
		return !p.isAlive();
	}

	/**
	 * Returns the result of the executed process. MUST only be called if {@link #isFinished()} returns true.
	 */
	public int getResultCode() {
		if (!isFinished()) {
			throw new IllegalStateException("Cannot get process result before it finishes!");
		}

		return p.exitValue();
	}

	public void waitFor() throws InterruptedException {
		p.waitFor();
	}
}
