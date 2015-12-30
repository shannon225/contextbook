package edu.washington.gs.maccoss.encyclopedia.utils.threading;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

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
		StringBuilder sb=new StringBuilder();
		OS os=OSDetector.getOS();
		if (OS.WINDOWS==os) {
			sb.append("cmd /c");
		}
		
		for (String arg : cmdAndArgs) {
			if (sb.length()>0) {
				sb.append(" ");
			}
			sb.append(arg);
		}
		Logger.logLine("Executing ["+sb.toString()+"]");
		ProcessBuilder pb=new ProcessBuilder(cmdAndArgs);
		p=pb.start();

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
	
	public void waitFor() throws InterruptedException {
		p.waitFor();
	}
}
