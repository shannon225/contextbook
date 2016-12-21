package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.io.FilterInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.ConcurrentModificationException;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class ProgressInputStream extends FilterInputStream {
	private static final int UPDATE_THRESHOLD=65536;
	
	private final ArrayList<ChangeListener> listeners=new ArrayList<ChangeListener>();
	private volatile long mark=0;
	private volatile long lastTriggeredLocation=0;
	private volatile long location=0;

	public ProgressInputStream(InputStream in) {
		super(in);
	}

	public void addChangeListener(ChangeListener l) {
		if (!listeners.contains(l)) {
			listeners.add(l);
		}
	}

	public void removeChangeListener(ChangeListener l) {
		listeners.remove(l);
	}

	public long getProgress() {
		return location;
	}

	@Override
	public int read() throws IOException {
		int i=super.read();
		if (i!=-1) {
			triggerChanged(location++);
		}
		return i;
	}

	@Override
	public int read(byte[] b, int off, int len) throws IOException {
		int i=super.read(b, off, len);
		if (i>0) {
			triggerChanged(location+=i);
		}
		return i;
	}

	@Override
	public long skip(long n) throws IOException {
		long i=super.skip(n);
		if (i>0) {
			triggerChanged(location+=i);
		}
		return i;
	}

	@Override
	public void mark(int readlimit) {
		super.mark(readlimit);
		mark=location;
	}

	@Override
	public void reset() throws IOException {
		super.reset();
		if (location!=mark) {
			triggerChanged(location=mark);
		}
	}

	private void triggerChanged(long location) {
		if (UPDATE_THRESHOLD>location-lastTriggeredLocation) {
			return;
		}
		lastTriggeredLocation=location;
		if (listeners.size()<=0) {
			return;
		}
		try {
			ChangeEvent evt=new ChangeEvent(this);
			for (ChangeListener l : listeners) {
				l.stateChanged(evt);
			}
		} catch (ConcurrentModificationException e) {
			triggerChanged(location);
		}
	}
}
