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

import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneAuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizationScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.ScoringResultsToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class ChromatogramAlignedPhosphoSiteLocalizer {
	// TODO add saving / exporting of library results. This will have to go into the threaded task (and an object saver)
	// TODO add execution job that can be queued
	// TODO add actual writing of library results to tables
	// TODO write a main class
	
	public static ArrayList<PercolatorPeptide> analyzeLocalization(ProgressIndicator progress, final SearchJobData job, boolean limitToQuantifiable, ArrayList<PercolatorPeptide> localPassingPSMIDs, final Optional<PeakLocationInferrer> inferrer, StripeFileInterface stripeFile, LibraryInterface searchedLibrary, final SearchParameters parameters) {
		HashMap<String, PSMData> uniquedData=parseDetectedPeptides(job, localPassingPSMIDs, inferrer, parameters);
		
		try {
			PhosphoLocalizer localizer=new PhosphoLocalizer(stripeFile, searchedLibrary, parameters);
			BackgroundFrequencyInterface background=localizer.getBackground();
			
			int cores=parameters.getNumberOfThreadsUsed();
			Logger.logLine("Extracting "+uniquedData.size()+" peptides...");

			Logger.logLine("Processing precursors scans...");
			String filename=stripeFile.getOriginalFileName();	
			File featureFile=new File(stripeFile.getFile().getAbsolutePath()+".localization.feature.txt");
			File outputFile=new File(stripeFile.getFile().getAbsolutePath()+".localization.report.txt");
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
			PeptideScoringResultsConsumer resultsConsumer=new ScoringResultsToTSVConsumer(featureFile, stripeFile.getFile(), EncyclopediaOneAuxillaryPSMScorer.getScoreNames(false), resultsQueue, 1);
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
						AnnotatedLibraryEntry entry=FragmentationModel.generateEntry(psm.getPeptideModSeq(), filename, psm.getAccessions(), psm.getPrecursorCharge(), psm.getRetentionTime(), parameters);
						String smartDecoy=PeptideUtils.reverse(psm.getPeptideModSeq(), parameters);
						FastaPeptideEntry decoyPeptide=new FastaPeptideEntry(filename, LibraryEntry.DECOY_STRING+psm.getAccession(), smartDecoy);
						AnnotatedLibraryEntry decoy=FragmentationModel.generateEntry(smartDecoy, decoyPeptide.getFilename(), decoyPeptide.getAccessions(), psm.getPrecursorCharge(), psm.getRetentionTime(), parameters);
						
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
			ArrayList<PercolatorPeptide> passingPeptides=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorVersionNumber(), featureFile, outputFile, parameters.getEffectivePercolatorThreshold());
			stripeFile.close();
			
			return passingPeptides;
			
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
