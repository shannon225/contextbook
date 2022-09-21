package edu.washington.gs.maccoss.encyclopedia.jobs;

import java.util.ArrayList;

public interface JobProcessor {

	ArrayList<WorkerJob> getQueue();

	void addJob(WorkerJob job);

	void fireJobUpdated(WorkerJob job);

}