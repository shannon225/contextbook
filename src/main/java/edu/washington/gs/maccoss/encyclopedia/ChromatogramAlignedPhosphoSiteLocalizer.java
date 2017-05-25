package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Optional;
import java.util.Vector;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ParsimonyProteinGrouper;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.AlternatePeakLocationInferrer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneAuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizationJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizationScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.ScoringResultsToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableConcatenator;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class ChromatogramAlignedPhosphoSiteLocalizer {
	// TODO add saving / exporting of library results. This will have to go into the threaded task (and an object saver)
	// TODO add execution job that can be queued
	// TODO add actual writing of library results to tables

	public static void convert(HashMap<String, String> arguments) {
		if (!arguments.containsKey("-i")||!arguments.containsKey("-l")||!arguments.containsKey("-o")) {
			Logger.errorLine("You are required to specify an input file or directory (-i), an input library file (-l) and an output library file (-o)");
			System.exit(1);
		}

		File diaFile=new File(arguments.get("-i"));
		File libraryFile=new File(arguments.get("-l"));
		File outputFile=new File(arguments.get("-o"));
		boolean alignBetweenFiles=SearchParameterParser.getBoolean("-a", arguments, true);
		boolean writeBlib=SearchParameterParser.getBoolean("-blib", arguments, false);
		
		SearchParameters parameters=SearchParameterParser.parseParameters(arguments);
		LibraryScoringFactory factory=new EncyclopediaOneScoringFactory(parameters);
		Logger.timelessLogLine("SearchToLIB EncyclopeDIA version "+factory.getVersion());

		Logger.timelessLogLine("Parameters:");
		Logger.timelessLogLine(" -i "+diaFile.getAbsolutePath());
		Logger.timelessLogLine(" -l "+libraryFile.getAbsolutePath());
		Logger.timelessLogLine(" -o "+outputFile.getAbsolutePath());
		Logger.timelessLogLine(" -a "+alignBetweenFiles);
		Logger.timelessLogLine(" -blib "+writeBlib);
		Logger.timelessLogLine(parameters.toString());
		
		EmptyProgressIndicator progress=new EmptyProgressIndicator();
		try {
			LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
			
			ArrayList<SearchJobData> pecanJobs=new ArrayList<SearchJobData>();
			if (diaFile.isDirectory()) {
				File[] files=diaFile.listFiles(StripeFileGenerator.getFilenameFilter());
				if (files.length==0) {
					Logger.errorLine("Your specified input (-i) directory didn't contain any .RAW files!");
					System.exit(1);
				}
				for (File file : files) {
					EncyclopediaJobData job=new EncyclopediaJobData(file, library, factory);
					pecanJobs.add(job);
				}
			} else {
				EncyclopediaJobData job=new EncyclopediaJobData(diaFile, library, factory);
				pecanJobs.add(job);
			}
			convert(progress, pecanJobs, outputFile, alignBetweenFiles, writeBlib, parameters);
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}

	public static void convert(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, File libFile, boolean alignBetweenFiles, boolean writeBlib, SearchParameters parameters) {
		Logger.logLine("Attempting to localize "+pecanJobs.size()+" searches...");
		ArrayList<PhosphoLocalizationJobData> phosphoJobs=new ArrayList<PhosphoLocalizationJobData>();
		for (SearchJobData job : pecanJobs) {
			PhosphoLocalizationJobData phosphoJob=analyzeLocalization(progress, job, false, Optional.empty(), parameters);
			phosphoJobs.add(phosphoJob);
		}
		
		Logger.logLine("Attempting to write "+pecanJobs.size()+" searches...");
		convert(progress, phosphoJobs, libFile, writeBlib, alignBetweenFiles);
	}
	
	public static void convert(ProgressIndicator progress, ArrayList<PhosphoLocalizationJobData> originalJobList, File libFile, boolean writeBlib, boolean alignBetweenFiles) {
		ArrayList<SearchJobData> processedJobs=new ArrayList<SearchJobData>();
		ArrayList<File> featureFiles=new ArrayList<File>();
		SearchJobData representativeJob=null;
		for (int i=0; i<originalJobList.size(); i++) {
			SearchJobData job=originalJobList.get(i);
			if (!job.hasBeenRun()) {
				Logger.logLine("Can't find a "+job.getSearchType()+" analysis of "+job.getDiaFile().getName()+", skipping extraction on that file.");
				continue;
			} else {
				processedJobs.add(job);
			}
			if (representativeJob==null) {
				representativeJob=job;
			}
			featureFiles.add(job.getFeatureFile());
		}

		if (representativeJob==null) {
			Logger.errorLine("Can't find any representative jobs! Failing...");

			for (int i=0; i<processedJobs.size(); i++) {
				SearchJobData job=processedJobs.get(i);
				Logger.errorLine(" Checking raw file "+(i+1)+": "+job.getDiaFile().exists());
				Logger.errorLine(" Checking feature file "+(i+1)+": "+job.getFeatureFile().exists());
				Logger.errorLine(" Checking result file "+(i+1)+": "+job.getOutputFile().exists());
			}
			return;
		}
		Logger.logLine("Using "+representativeJob.getDiaFile().getName()+" to extract representative search parameters");

		String filename=libFile.getName();
		if (filename.lastIndexOf('.')>0) {
			filename=filename.substring(0, filename.lastIndexOf('.'));
		}
		File bigFeatureFile=new File(representativeJob.getFeatureFile().getParentFile(), filename+"_concatenated_features.txt");
		File bigPercolatorFile=new File(representativeJob.getFeatureFile().getParentFile(), filename+"_concatenated_results.txt");
		
		SearchParameters parameters=representativeJob.getParameters();
		float threshold=parameters.getEffectivePercolatorThreshold();
		try {
			ArrayList<PercolatorPeptide> passingPeptides;
			if (featureFiles.size()==1) {
				// if there's only one file then don't need to re-run percolator
				passingPeptides=PercolatorReader.getPassingPeptidesFromTSV(representativeJob.getOutputFile(), threshold);
			} else if (bigPercolatorFile.exists()&&bigPercolatorFile.canRead()) {
				// if we've already run percolator then don't need to re-run percolator
				passingPeptides=PercolatorReader.getPassingPeptidesFromTSV(bigPercolatorFile, threshold);
			} else {
				TableConcatenator.concatenateTables(featureFiles, bigFeatureFile);
				passingPeptides=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorVersionNumber(), bigFeatureFile, bigPercolatorFile, threshold);
			}
			
			ArrayList<ProteinGroup> proteins=ParsimonyProteinGrouper.groupProteins(passingPeptides);
			Logger.logLine("Identified "+passingPeptides.size()+" peptides ("+proteins.size()+" proteins) across all files at a "+(threshold*100.0f)+"% FDR threshold.");

			Optional<PeakLocationInferrer> inferrer;
			if (alignBetweenFiles) {
				if (processedJobs.size()>1) {
					Logger.logLine("Inferring peak boundaries across files...");
					inferrer=Optional.of(AlternatePeakLocationInferrer.getAlignmentData(new EmptyProgressIndicator(), processedJobs, passingPeptides, parameters));
					Logger.logLine("...Finished peak inference.");
				} else {
					Logger.logLine("Only processing one file so no peak inference is necessary.");
					inferrer=Optional.empty();
				}
			} else {
				Logger.logLine("User requested no RT alignment between files.");
				inferrer=Optional.empty();
			}

			if (writeBlib) {
				SearchToBLIB.convertBlib(progress, processedJobs, libFile, Optional.of(passingPeptides), inferrer);
			} else {
				SearchToBLIB.convertElib(progress, processedJobs, libFile, Optional.of(passingPeptides), inferrer, proteins, parameters);
			}
			
			progress.update(passingPeptides.size()+" peptides identified at "+(threshold*100.0f)+"% FDR", 1.0f);
		} catch (IOException ioe) {
			Logger.errorLine("Error creating concatenated feature file");
			Logger.errorException(ioe);
		} catch (InterruptedException ie) {
			Logger.errorLine("Error creating concatenated feature file");
			Logger.errorException(ie);
		}
	}
	
	public static PhosphoLocalizationJobData analyzeLocalization(ProgressIndicator progress, final SearchJobData job, boolean limitToQuantifiable, final Optional<PeakLocationInferrer> inferrer, final SearchParameters parameters) {
		StripeFileInterface stripeFile=StripeFileGenerator.getFile(job.getDiaFile(), parameters);		
		ArrayList<PercolatorPeptide> localPassingPSMIDs=PercolatorReader.getPassingPeptidesFromTSV(job.getOutputFile(), parameters.getEffectivePercolatorThreshold());

		HashMap<String, String> paramMap=parameters.toParameterMap();
		paramMap.put("-runPhosphoLocalization", "false");
		PhosphoLocalizationJobData phosphoJob=new PhosphoLocalizationJobData(stripeFile.getFile(), SearchParameterParser.parseParameters(paramMap));
		String filename=stripeFile.getOriginalFileName();	
		File featureFile=phosphoJob.getFeatureFile();
		File outputFile=phosphoJob.getOutputFile();
		
		HashMap<String, PSMData> uniquedData=parseDetectedPeptides(job, localPassingPSMIDs, inferrer, parameters);
		
		try {
			PhosphoLocalizer localizer=new PhosphoLocalizer(stripeFile, parameters);
			BackgroundFrequencyInterface background=localizer.getBackground();
			
			int cores=parameters.getNumberOfThreadsUsed();
			Logger.logLine("Localizing "+uniquedData.size()+" peptides...");

			Logger.logLine("Processing precursors scans...");
			PrecursorScanMap precursors=new PrecursorScanMap(stripeFile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));
			
			// get targeted ranges
			ArrayList<Range> ranges=new ArrayList<Range>();
			for (Range range : stripeFile.getRanges().keySet()) {
				if (!parameters.useTargetWindowCenter()||range.contains(parameters.getTargetWindowCenter())) {
					ranges.add(range);
				}
			}
			Collections.sort(ranges);

			BlockingQueue<PeptideScoringResult> resultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
			PeptideScoringResultsConsumer resultsConsumer=new ScoringResultsToTSVConsumer(featureFile, stripeFile.getFile(), EncyclopediaOneAuxillaryPSMScorer.getScoreNames(true), resultsQueue, 1);
			Thread consumerThread=new Thread(resultsConsumer);
			consumerThread.start();
			
			// get stripes
			int rangesFinished=0;
			float numberOfTasks=2.0f+ranges.size();
			float baseIncrement=1.0f/numberOfTasks;
			
			for (Range range : ranges) {
				String baseMessage="Working on "+range+" m/z";
				float baseProgress=(1.0f+rangesFinished)/numberOfTasks;
				progress.update(baseMessage, baseProgress);

				Logger.logLine("Processing "+range);
				
				LibraryBackgroundInterface libraryBackground=background.getLibraryBackground(range.getMiddle(), parameters.getFragmentTolerance());
				EncyclopediaScorer scorer=new EncyclopediaOneScorer(parameters, libraryBackground);

				ArrayList<Stripe> stripes=stripeFile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, false);
				Collections.sort(stripes);

				// prepare executor for background
				ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("STRIPE_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
				LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
				ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

				int count=0;
				for (PSMData psm : uniquedData.values()) {
					if (range.contains((float)psm.getPrecursorMZ())) {
						count++;
						AnnotatedLibraryEntry entry=FragmentationModel.generateEntry(psm.getPeptideModSeq(), filename, psm.getAccessions(), psm.getPrecursorCharge(), psm.getRetentionTime(), false, parameters);
						String smartDecoy=PeptideUtils.reverse(psm.getPeptideModSeq(), parameters);
						FastaPeptideEntry decoyPeptide=new FastaPeptideEntry(filename, LibraryEntry.DECOY_STRING+psm.getAccession(), smartDecoy);
						AnnotatedLibraryEntry decoy=FragmentationModel.generateEntry(smartDecoy, decoyPeptide.getFilename(), decoyPeptide.getAccessions(), psm.getPrecursorCharge(), psm.getRetentionTime(), true, parameters);
						
						ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
						entries.add(entry);
						entries.add(decoy);
						
						executor.submit(new PhosphoLocalizationScoringTask(scorer, entries, stripes, precursors, localizer, resultsQueue, parameters));
					}
				}

				executor.shutdown();
				while (!executor.isTerminated()) {
					Logger.logLine(workQueue.size()+" peptides remaining for "+range+"...");
					float finishedFraction=(count-workQueue.size())/(float)count;
					progress.update(baseMessage, baseProgress+baseIncrement*(0.2f+finishedFraction*0.8f));
					Thread.sleep(500);
				}
				executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
				
				rangesFinished++;
			}

			resultsQueue.put(PeptideScoringResult.POISON_RESULT);

			consumerThread.join();
			resultsConsumer.close();
			Logger.logLine("Finished generating feature file, analyzed "+resultsConsumer.getNumberProcessed()+" peptides.");

			progress.update("Running Percolator", (1.0f+rangesFinished)/numberOfTasks);
			PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorVersionNumber(), featureFile, outputFile, parameters.getEffectivePercolatorThreshold());
			stripeFile.close();
			
			return phosphoJob;
			
		} catch (IOException ioe) {
			Logger.errorLine("Error processing "+stripeFile.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ioe);
		} catch (SQLException sqle) {
			Logger.errorLine("Error processing "+stripeFile.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", sqle);
		} catch (DataFormatException dfe) {
			Logger.errorLine("Error processing "+stripeFile.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", dfe);
		} catch (InterruptedException ie) {
			Logger.errorLine("Error processing "+stripeFile.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ie);
		}
	}

	private static HashMap<String, PSMData> parseDetectedPeptides(final SearchJobData job, ArrayList<PercolatorPeptide> localPassingPSMIDs, final Optional<PeakLocationInferrer> inferrer,
			final SearchParameters parameters) {
		Logger.logLine("Number of peptides: "+localPassingPSMIDs.size());

		final TObjectFloatHashMap<String> localSavedIDs=new TObjectFloatHashMap<String>();
		for (PercolatorPeptide psm : localPassingPSMIDs) {
			localSavedIDs.put(psm.getPsmID(), psm.getQValue());
		}

		final Vector<PSMData> data=new Vector<PSMData>();
		
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String psmID=row.get("id");
				String peptideModSeq=PercolatorPeptide.getPeptideSequence(psmID);
					boolean isDecoy=PercolatorPeptide.isPSMIDDecoy(psmID);
					if (!isDecoy) {
						float retentionTime;// in seconds
						int scanID;
						
						// prefer actual identification, fall back on RT inference
						if (localSavedIDs.contains(psmID)) {
							String rtString=row.get("midTime"); // in seconds
							if (rtString!=null) {
								retentionTime=Float.parseFloat(rtString);
							} else {
								rtString=row.get("RTinMin"); // in minutes so *60
								retentionTime=Float.parseFloat(rtString)*60f;
							}
							scanID=Integer.parseInt(row.get("ScanNr"));
						} else if (inferrer.isPresent()) {
							if (localSavedIDs.contains(psmID)) {
								String rtString=row.get("midTime"); // in seconds
								if (rtString!=null) {
									retentionTime=Float.parseFloat(rtString);
								} else {
									rtString=row.get("RTinMin"); // in minutes so *60
									retentionTime=Float.parseFloat(rtString)*60f;
								}
								
								float warpedRT=inferrer.get().getPreciseRTInSec(job, peptideModSeq, retentionTime);
								if (warpedRT!=retentionTime) {
									Logger.errorLine("Don't trust ID for "+peptideModSeq+" (global RT:"+warpedRT+", local RT:"+retentionTime+"). Using the warped RT!");
									// warping is better (original is way outside the warping margins)
									retentionTime=warpedRT;
									scanID=-1; // negative scan ID for inferred IDs
								} else {
									// original detection is better (within the warping margins)
									scanID=Integer.parseInt(row.get("ScanNr"));
								}
									
							} else {
								// no detect, so use warped retention time
								retentionTime=inferrer.get().getWarpedRTInSec(job, peptideModSeq);
								scanID=-1; // negative scan ID for inferred IDs
							}
						} else {
							if (localSavedIDs.contains(psmID)) {
								// no warping, so use RT
								String rtString=row.get("midTime"); // in seconds
								if (rtString!=null) {
									retentionTime=Float.parseFloat(rtString);
								} else {
									rtString=row.get("RTinMin"); // in minutes so *60
									retentionTime=Float.parseFloat(rtString)*60f;
								}
								scanID=Integer.parseInt(row.get("ScanNr"));
							} else {
								// not in local search and no warping available
								return;
							}
						}

						double precursorMZ=Double.parseDouble(row.get("precursorMz"));
						// FIXME need to get peptide charge from window
						byte precursorCharge=PercolatorPeptide.getCharge(psmID);

						float sortingScore;
						String sortingScoreString=row.get("primary"); // Encyclopedia/XCordia
						if (sortingScoreString==null) {
							sortingScoreString=row.get("xTandem"); // old Encyclopedia
						}
						if (sortingScoreString==null) {
							sortingScoreString=row.get("peakZScore"); // Pecan
						}
						if (sortingScoreString==null) {
							sortingScoreString=row.get("peakBGScore"); // Pecan
						}
						if (sortingScoreString==null) {
							Logger.errorLine("Can't parse score from header from ["+row.keySet()+"]");
							throw new EncyclopediaException("Can't parse score from header from ["+row.keySet()+"]");
						}
						sortingScore=Float.parseFloat(sortingScoreString);

						String samplingTimeString=row.get("sampledTimes");
						float duration=samplingTimeString==null?(parameters.getExpectedPeakWidth()):Float.parseFloat(samplingTimeString);

						String proteinString=row.get("protein");
						HashSet<String> accessions=PSMData.stringToAccessions(proteinString);
						data.add(new PSMData(accessions, scanID, precursorMZ, precursorCharge, peptideModSeq, retentionTime, sortingScore, sortingScore, duration));
					}
				}
		};
		
		TableParser.parseTSV(job.getFeatureFile(), muscle);

		Logger.logLine("Parsed features and scores for "+data.size()+" peptides.");
		HashMap<String, PSMData> uniquedData=new HashMap<String, PSMData>();
		for (PSMData psmData : data) {
			String key=psmData.getPeptideModSeq()+"+"+psmData.getPrecursorCharge();
			PSMData prev=uniquedData.get(key);
			if (prev!=null) {
				if (prev.getSortingScore()<psmData.getSortingScore()) {
					// scores scores are high
					uniquedData.put(key, psmData);
				}
			} else {
				uniquedData.put(key, psmData);
			}
		}
		return uniquedData;
	}
}
