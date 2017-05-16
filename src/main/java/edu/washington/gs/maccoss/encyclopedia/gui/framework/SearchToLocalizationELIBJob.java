package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.ChromatogramAlignedPhosphoSiteLocalizer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;

public class SearchToLocalizationELIBJob extends SwingJob {
	private final File elibFile;
	private final boolean alignBetweenFiles;
	private final SearchParameters parameters;

	public SearchToLocalizationELIBJob(File blibFile, boolean alignBetweenFiles, SearchParameters parameters, JobProcessor processor) {
		super(processor);
		this.elibFile=blibFile;
		this.alignBetweenFiles=alignBetweenFiles;
		this.parameters=parameters;
	}
	
	@Override
	public String getJobTitle() {
		return "Phospho localizing "+elibFile.getName();
	}

	@Override
	public void runJob() throws Exception {
		ArrayList<SearchJobData> jobData=new ArrayList<SearchJobData>();
		for (SwingJob job : processor.getQueue()) {
			if (job instanceof SearchJob) {
				jobData.add(((SearchJob)job).getSearchData());
			}
		}

		ChromatogramAlignedPhosphoSiteLocalizer.convert(getProgressIndicator(), jobData, elibFile, false, alignBetweenFiles, parameters);
	}
}
