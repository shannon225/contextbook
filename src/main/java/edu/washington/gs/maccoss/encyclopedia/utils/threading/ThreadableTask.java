package edu.washington.gs.maccoss.encyclopedia.utils.threading;

import java.util.concurrent.Callable;

public abstract class ThreadableTask<R> implements Callable<R> {
	public abstract String getTaskName();
	
	protected abstract R process();
	
	public ThreadableTask() {
	}
	
	@Override
	public R call() {
		//Logger.logLine("Started "+getTaskName());
		//long startTime=System.currentTimeMillis();

		R ret=process();

		//long stopTime=System.currentTimeMillis();
		//Logger.logLine("Finished "+getTaskName()+" in "+(stopTime-startTime)/1000+" seconds.");
		return ret;
	}
}
