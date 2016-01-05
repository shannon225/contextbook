package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.util.ArrayList;

public interface JobProcessor {

	ArrayList<SwingJob> getQueue();

	void addJob(SwingJob job);

	void fireJobUpdated(SwingJob job);

}