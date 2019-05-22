package edu.washington.gs.maccoss.encyclopedia.utils.threading;

import java.util.concurrent.LinkedBlockingQueue;

/**
 * basically just turns offer (submit job, ignore result) into put (submit job and wait until there's space) unless interrupted. 
 *
 * @param <E>
 */
public class LimitedQueue<E> extends LinkedBlockingQueue<E> {
	private static final long serialVersionUID = 1L;

	public LimitedQueue(int capacity) {
		super(capacity);
	}

	@Override
	public boolean offer(E e) {
		try {
			put(e);
			return true;
		} catch (InterruptedException ie) {
			Thread.currentThread().interrupt();
		}
		return false;
	}
	
	
}
