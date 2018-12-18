package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import java.io.File;
import java.util.Optional;

import javax.swing.ImageIcon;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;

public interface ParametersPanelInterface {
	ProgramType getProgram();
	ImageIcon getImage();
	
	Optional<String> canLoadData();
	SwingJob getJob(File diaFile, JobProcessorTableModel model);
	void askForSetupFile();
	SearchParameters getParameters();
	void savePreferences();
	File getBackgroundFastaFile();
}