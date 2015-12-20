package edu.washington.gs.maccoss.encyclopedia.utils.threading;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class ExternalExecutor {
	private final String[] cmdAndArgs;
	private volatile Process p;
	
	public ExternalExecutor(String[] cmdAndArgs) {
		this.cmdAndArgs=cmdAndArgs;
	}

	public BlockingQueue<String> start() throws IOException {
		final BlockingQueue<String> queue=new LinkedBlockingQueue<String>();
		StringBuilder sb=new StringBuilder();
		for (String arg : cmdAndArgs) {
			if (sb.length()>0) {
				sb.append(" ");
			}
			sb.append(arg);
		}
		Logger.logLine("Executing ["+sb.toString()+"]");
		p=Runtime.getRuntime().exec(cmdAndArgs);

		new Thread(new Runnable() {
			public void run() {
				BufferedReader input=new BufferedReader(new InputStreamReader(p.getInputStream()));
				String line=null;

				try {
					while ((line=input.readLine())!=null) {
						System.out.println("--> "+line);
						queue.add(line);
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
	
	public void waitFor() throws InterruptedException {
		p.waitFor();
	}
}
