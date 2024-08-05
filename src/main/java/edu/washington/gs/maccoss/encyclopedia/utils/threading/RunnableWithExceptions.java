package edu.washington.gs.maccoss.encyclopedia.utils.threading;

import java.util.Optional;

public abstract class RunnableWithExceptions implements Runnable {
	Exception e=null;
	
	protected void setException(Exception e) {
		this.e=e;
	}

	public Optional<Exception> getException() {
		return Optional.ofNullable(e);
	}
}
