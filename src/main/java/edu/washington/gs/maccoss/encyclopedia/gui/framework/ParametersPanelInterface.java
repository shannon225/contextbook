package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import java.io.File;
import java.util.Optional;

import javax.swing.ImageIcon;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;

public interface ParametersPanelInterface {
	String getProgramName();
	String getCitation();
	String getAboutMessage();
	ImageIcon getImage();
	
	Optional<String> canLoadData();
	SwingJob getJob(File diaFile, JobProcessorTableModel model);
	void askForSetupFile();
	SearchParameters getParameters();
	void savePreferences();
}