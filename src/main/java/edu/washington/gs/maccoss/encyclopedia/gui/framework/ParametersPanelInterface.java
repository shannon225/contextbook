package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import java.io.File;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;

public interface ParametersPanelInterface {

	Optional<String> canLoadData();
	SwingJob getJob(File diaFile, JobProcessorTableModel model);

}