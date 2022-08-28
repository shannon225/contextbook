package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.SearchToBLIB;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.EncyclopediaTwoPeakLocationInferrer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.DIAProcessor;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.spectrumprocessors.WindowDownsampler;
import edu.washington.gs.maccoss.encyclopedia.filewriters.LibraryUtilities;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class CombineELIBsAndExtractGroupSpecificLibrariesJob extends SwingJob {
	private final File saveDirectory;
	private final Optional<File> singleInjectionExample;
	private final SearchParameters parameters;

	public CombineELIBsAndExtractGroupSpecificLibrariesJob(File saveDirectory, Optional<File> singleInjectionExample, SearchParameters parameters, JobProcessor processor) {
		super(processor);
		this.saveDirectory=saveDirectory;
		this.singleInjectionExample=singleInjectionExample;
		this.parameters=parameters;
	}
	
	@Override
	public String getJobTitle() {
		return "Write Batch-specific Libraries to "+saveDirectory.getName();
	}

	@Override
	public void runJob() throws Exception {
		ProgressIndicator progress = getProgressIndicator();
		// make sure saveDirectory is indeed a directory
		if (saveDirectory.exists()) {
			if(!saveDirectory.isDirectory()) {
				Logger.errorLine("You need to specify a save directory (not a file)!");
				return;
			}
		} else {
			if (!saveDirectory.mkdirs()) {
				Logger.errorLine("Problem creating directory ["+saveDirectory.getAbsolutePath()+"]!");
				return;
			}
		}
		
		WindowDownsampler downsampler=null;
		if (singleInjectionExample.isPresent()) {
			StripeFileInterface stripeFile=StripeFileGenerator.getFile(singleInjectionExample.get(), parameters, true);
			ArrayList<Range> downsampledRanges=new ArrayList<Range>(stripeFile.getRanges().keySet());
			Collections.sort(downsampledRanges);
			stripeFile.close();
			Logger.logLine("Downsampling quant data to "+downsampledRanges.size()+" windows");

			downsampler=new WindowDownsampler(downsampledRanges, parameters.getFragmentTolerance());
		}
		
		// grab jobs from the current queue and downsample DIA data
		
		ArrayList<SearchJobData> jobData=new ArrayList<SearchJobData>();
		ArrayList<SwingJob> queue = processor.getQueue();
		for (SwingJob job : queue) {
			if (job instanceof SearchJob) {
				SearchJobData searchData = ((SearchJob)job).getSearchData();
				
				if (downsampler!=null) {
					DIAProcessor processor=new DIAProcessor(downsampler, parameters);
					File originalFile=searchData.getDiaFileReader().getFile();
					File newQuantFile=new File(originalFile.getParentFile(), originalFile.getName()+".downsampled"+StripeFile.DIA_EXTENSION);
					Logger.logLine("Downsampling "+originalFile.getName()+" to create "+newQuantFile.getName());
					processor.processStripeFile(new SubProgressIndicator(progress, 0.25f/queue.size()), searchData.getDiaFileReader(), newQuantFile, false);
					searchData=searchData.updateQuantFile(newQuantFile);
				}
				jobData.add(searchData);
			}
		}
		Logger.logLine("Found "+jobData.size()+" jobs in the queue to combine...");
		
		if (jobData.size()<=1) {
			Logger.errorLine("You need to queue more than one batch-specific libraries to combine them!");
			return;
		}
		
		// force the max number of quant peaks to be all peaks!
		HashMap<String, String> params=parameters.toParameterMap();
		params.put(SearchParameters.NUMBER_OF_QUANTITATIVE_PEAKS, Integer.toString(Integer.MAX_VALUE));
		SearchParameters quantParameters=SearchParameterParser.parseParameters(params);
		
		// combine files
		progress.update("Calculating global FDR across batch-specific libraries");
		Logger.logLine("Calculating global FDR across batch-specific libraries");
		File intermediateQuantLibraryFile=new File(saveDirectory, "batch_combined_quant_report.elib");
		SearchToBLIB.convert(new SubProgressIndicator(progress, 0.25f), jobData, intermediateQuantLibraryFile, false, true, quantParameters);

		// identify bestQuant ions
		progress.update("Calculating global transitions for quantification");
		Logger.logLine("Calculating global transitions for quantification");
		Pair<HashMap<SearchJobData, TObjectFloatHashMap<String>>, HashMap<String, double[]>> archetypalData=EncyclopediaTwoPeakLocationInferrer.getArchetypals(new SubProgressIndicator(progress, 0.25f), jobData, parameters);
		HashMap<String, double[]> globalTransitions=archetypalData.y;
		int numWithAtLeastX=0;
		int numWithAtLeastY=0;
		int max=0;
		for (double[] targets : globalTransitions.values()) {
			if (targets.length>=parameters.getMinNumOfQuantitativePeaks()) numWithAtLeastX++;
			if (targets.length>=parameters.getNumberOfQuantitativePeaks()) numWithAtLeastY++;
			if (targets.length>max) max=targets.length;
		}
		Logger.logLine("Found "+globalTransitions.size()+"/"+numWithAtLeastX+"/"+numWithAtLeastY+" peptides with at least 0/"+parameters.getMinNumOfQuantitativePeaks()+"/"+parameters.getNumberOfQuantitativePeaks()+" targetable fragment ions");
		Logger.logLine("Max number of transitions considered: "+max);
		
		// then part them out
		progress.update("Extract individual batch libraries from global analysis");
		Logger.logLine("Extract individual batch libraries from global analysis");
		LibraryFile library=new LibraryFile();
		library.openFile(intermediateQuantLibraryFile);
		LibraryUtilities.extractSampleSpecificLibraries(new SubProgressIndicator(progress, 0.25f), saveDirectory, Optional.of(globalTransitions), library, parameters);
	}
}
