package edu.washington.gs.maccoss.encyclopedia.gui.pecan;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.SwingWorker;

import com.google.common.base.Optional;

import edu.washington.gs.maccoss.encyclopedia.Pecanpie;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.io.ProgressMessage;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class PecanJob extends SwingWorker<Nothing, ProgressMessage> {
	private final Optional<ArrayList<FastaEntry>> targetList;
	private final File diaFile;
	private final File fastaFile;
	private final File featureFile;
	private final File outputFile;
	private final PecanScoringFactory taskFactory;

	private volatile String message="";
	private volatile float progress=0.0f;
	private final PecanFileProcessorModel processor;

	public PecanJob(PecanFileProcessorModel processor, Optional<ArrayList<FastaEntry>> targetList, File diaFile, File fastaFile, File featureFile, File outputFile, PecanScoringFactory taskFactory) {
		this.targetList=targetList;
		this.diaFile=diaFile;
		this.fastaFile=fastaFile;
		this.featureFile=featureFile;
		this.outputFile=outputFile;
		this.taskFactory=taskFactory;
		this.processor=processor;
	}

	@Override
	protected void done() {
		message="Done!";
		progress=1.0f;
		processor.fireJobUpdated(this);
	}

	@Override
	protected Nothing doInBackground() throws Exception {
		try {
			final ProgressIndicator indicator=new ProgressIndicator() {
				@Override
				public void update(String message, float totalProgress) {
					publish(new ProgressMessage(message, totalProgress));
				}
			};
			
			Pecanpie.runPie(indicator, targetList, diaFile, fastaFile, featureFile, outputFile, taskFactory);
		} catch (Exception e) {
			publish(new ProgressMessage("Encountered Fatal Error!", -1.0f));
			e.printStackTrace();
		}
		return Nothing.NOTHING;
	}

	@Override
	protected void process(List<ProgressMessage> chunks) {
		for (ProgressMessage p : chunks) {
			if (progress<=p.getProgress()) {
				progress=p.getProgress();
				message=p.getMessage();
				processor.fireJobUpdated(this);
			}
		}
	}
	
	public File getDiaFile() {
		return diaFile;
	}
	
	public File getOutputFile() {
		return outputFile;
	}
	
	public File getFeatureFile() {
		return featureFile;
	}
	
	public PecanScoringFactory getTaskFactory() {
		return taskFactory;
	}
	
	public ProgressMessage getProgressMessage() {
		return new ProgressMessage(message, progress);
	}
	
	public float getTotalProgress() {
		return progress;
	}
	
	public String getMessage() {
		return message;
	}
}
