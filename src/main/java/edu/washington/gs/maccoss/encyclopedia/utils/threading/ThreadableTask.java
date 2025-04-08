package edu.washington.gs.maccoss.encyclopedia.utils.threading;

import java.util.concurrent.Callable;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public abstract class ThreadableTask<R> implements Callable<R> {
	public abstract String getTaskName();
	
	protected abstract R process();
	
	public ThreadableTask() {
	}
	
	@Override
	public R call() {

		try {
			//Logger.logLine("Started "+getTaskName());
			//long startTime=System.currentTimeMillis();
			
			R ret=process();
			return ret;
			
			//long stopTime=System.currentTimeMillis();
			//Logger.logLine("Finished "+getTaskName()+" in "+(stopTime-startTime)/1000+" seconds.");
		} catch (Throwable t) {
			Logger.errorLine("Encountered unexpected exception!");
			t.printStackTrace();
			Logger.errorException(t);
			throw new EncyclopediaException("Encountered unexpected exception!", t);
		}

	}
}
