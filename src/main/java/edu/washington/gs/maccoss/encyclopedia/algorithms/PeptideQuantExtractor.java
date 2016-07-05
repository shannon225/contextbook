package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Optional;
import java.util.Vector;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class PeptideQuantExtractor {
	public static ArrayList<IntegratedLibraryEntry> parseSearchFeatures(ProgressIndicator progress, File f, boolean limitToQuantifiable, ArrayList<PercolatorPeptide> globalPassingPSMIDs, ArrayList<PercolatorPeptide> localPassingPSMIDs, StripeFileInterface stripeFile, LibraryInterface searchedLibrary, final SearchParameters parameters) {
		HashSet<String> passingPeptideSequences=new HashSet<String>();
		for (PercolatorPeptide psm : globalPassingPSMIDs) {
			String peptideModSeq=PercolatorPeptide.getPeptideSequence(psm.getPsmID());
			passingPeptideSequences.add(peptideModSeq);
		}
		
		final TObjectFloatHashMap<String> savedIDs=new TObjectFloatHashMap<String>();
		for (PercolatorPeptide psm : localPassingPSMIDs) {
			String peptideModSeq=PercolatorPeptide.getPeptideSequence(psm.getPsmID());
			if (passingPeptideSequences.contains(peptideModSeq)) {
				savedIDs.put(psm.getPsmID(), psm.getQValue());
			}
		}

		final Vector<PSMData> data=new Vector<PSMData>();
		
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String psmID=row.get("id");
				if (savedIDs.contains(psmID)) {
					boolean isDecoy=PercolatorPeptide.isPSMIDDecoy(psmID);
					if (!isDecoy) {
						
						int scanID=Integer.parseInt(row.get("ScanNr"));
						double precursorMZ=Double.parseDouble(row.get("precursorMz"));
						// FIXME need to get peptide charge from window
						byte precursorCharge=PercolatorPeptide.getCharge(psmID);
						String peptideModSeq=PercolatorPeptide.getPeptideSequence(psmID);

						float retentionTime;// in seconds
						String rtString=row.get("midTime"); // in seconds
						if (rtString!=null) {
							retentionTime=Float.parseFloat(rtString);
						} else {
							rtString=row.get("RTinMin"); // in minutes so *60
							retentionTime=Float.parseFloat(rtString)*60f;
						}
						float score=savedIDs.get(psmID);

						float sortingScore;
						String sortingScoreString=row.get("xTandem"); // Encyclopedia
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
						data.add(new PSMData(accessions, scanID, precursorMZ, precursorCharge, peptideModSeq, retentionTime, score, sortingScore, duration));
					}
				}
			}
		};
		
		TableParser.parseTSV(f, muscle);

		try {
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
			
			return extractPeptides(progress, searchedLibrary, stripeFile, uniquedData.values(), limitToQuantifiable, parameters);
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
	
	public static ArrayList<IntegratedLibraryEntry> extractPeptides(ProgressIndicator progress, LibraryInterface searchedLibrary, StripeFileInterface stripefile, Collection<PSMData> data, boolean limitToQuantifiable, SearchParameters parameters) throws IOException, SQLException, DataFormatException, InterruptedException {
		ConcurrentLinkedQueue<IntegratedLibraryEntry> savedEntries=new ConcurrentLinkedQueue<IntegratedLibraryEntry>();
		int cores=parameters.getNumberOfThreadsUsed();

		String filename=stripefile.getOriginalFileName();
		PhosphoLocalizer localizer=null;
		if (parameters.isRunPhosphoLocalization()&&searchedLibrary!=null) {
			localizer=new PhosphoLocalizer(stripefile, searchedLibrary, parameters);
		}
		
		// get identified peptides
		HashMap<String, PSMData> peptideModSeqs=new HashMap<String, PSMData>();
		for (PSMData psm : data) {
			peptideModSeqs.put(psm.getPeptideModSeq(), psm);
		}
		
		// get targeted ranges
		ArrayList<Range> ranges=new ArrayList<Range>();
		for (Range range : stripefile.getRanges().keySet()) {
			if (!parameters.useTargetWindowCenter()||range.contains(parameters.getTargetWindowCenter())) {
				ranges.add(range);
			}
		}
		Collections.sort(ranges);

		// get stripes
		int rangesFinished=0;
		float numberOfTasks=2.0f+ranges.size();
		for (Range range : ranges) {
			String baseMessage="Working on "+range+" m/z";
			float baseProgress=(1.0f+rangesFinished)/numberOfTasks;
			progress.update(baseMessage, baseProgress);

			Logger.logLine("Processing "+range);

			ArrayList<Stripe> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			Collections.sort(stripes);

			// prepare executor for background
			ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("SWATH_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
			LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
			ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 


			for (PSMData psm : data) {
				if (range.contains((float)psm.getPrecursorMZ())) {
					executor.submit(new PeptideQuantExtractorTask(filename, psm, Optional.ofNullable(localizer), stripes, parameters, savedEntries, limitToQuantifiable));
				}
			}

			executor.shutdown();
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
			
			rangesFinished++;
		}

		ArrayList<IntegratedLibraryEntry> entryList=new ArrayList<IntegratedLibraryEntry>();
		for (IntegratedLibraryEntry entry : savedEntries) {
			entryList.add(entry);
		}

		return entryList;
	}
}
