package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.jobs.WorkerJob;

public interface JobProcessor {

	ArrayList<WorkerJob> getQueue();

	void addJob(WorkerJob job);

	void fireJobUpdated(WorkerJob job);

}