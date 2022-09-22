package edu.washington.gs.maccoss.encyclopedia.jobs;

import java.util.ArrayList;

public class SimpleJobProcessor implements JobProcessor {
	private final ArrayList<WorkerJob> queue;
	

	@SuppressWarnings("unchecked")
	public SimpleJobProcessor(ArrayList<WorkerJob> queue) {
		this.queue = (ArrayList<WorkerJob>)queue.clone();
	}

	@SuppressWarnings("unchecked")
	@Override
	public ArrayList<WorkerJob> getQueue() {
		return (ArrayList<WorkerJob>)queue.clone();
	}

	@Override
	public void addJob(WorkerJob job) {
		queue.add(job);
	}

	@Override
	public void fireJobUpdated(WorkerJob job) {
		// ignore notifications
	}
	
}
