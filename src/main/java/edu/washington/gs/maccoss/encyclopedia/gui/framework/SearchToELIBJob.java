package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.SearchToBLIB;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;

public class SearchToELIBJob extends SwingJob {
	private final File elibFile;

	public SearchToELIBJob(File blibFile, JobProcessor processor) {
		super(processor);
		this.elibFile=blibFile;
	}
	
	@Override
	public String getJobTitle() {
		return "Write Library "+elibFile.getName();
	}

	@Override
	public void runJob() throws Exception {
		ArrayList<SearchJobData> jobData=new ArrayList<SearchJobData>();
		for (SwingJob job : processor.getQueue()) {
			if (job instanceof SearchJob) {
				jobData.add(((SearchJob)job).getSearchData());
			}
		}

		LibraryInterface libraryTemplate=null;
		SearchToBLIB.convert(getProgressIndicator(), jobData, elibFile, false, Optional.ofNullable(libraryTemplate));
	}
}
