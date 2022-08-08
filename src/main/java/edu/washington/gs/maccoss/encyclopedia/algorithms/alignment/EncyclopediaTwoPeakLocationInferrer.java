package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibrarySearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class EncyclopediaTwoPeakLocationInferrer {

	public static PeakLocationInferrerInterface getAlignmentData(ProgressIndicator progress, List<? extends SearchJobData> pecanJobs, ArrayList<PercolatorPeptide> passingPeptides, SearchParameters params) {
		ProgressIndicator subProgress1=new SubProgressIndicator(progress, 0.5f);
		Pair<HashMap<SearchJobData,TObjectFloatHashMap<String>>, HashMap<String,double[]>> pair=getArchetypals(subProgress1, pecanJobs, passingPeptides, params);
		return AlternatePeakLocationInferrer.getInferrer(progress, pecanJobs, pair, params);
	}

	public static Pair<HashMap<SearchJobData, TObjectFloatHashMap<String>>, HashMap<String, double[]>> getArchetypals(ProgressIndicator progress, List<? extends SearchJobData> jobs, SearchParameters params) {
		return getArchetypals(progress, jobs, null, params);
	}
	private static Pair<HashMap<SearchJobData, TObjectFloatHashMap<String>>, HashMap<String, double[]>> getArchetypals(ProgressIndicator progress, List<? extends SearchJobData> jobs, ArrayList<PercolatorPeptide> passingPeptides, SearchParameters params) {
		HashSet<String> passingPeptideModSeqs=new HashSet<>();
		// always accept if no passing peptides were used
		if (passingPeptides!=null) {
			for (PercolatorPeptide peptide : passingPeptides) {
				passingPeptideModSeqs.add(peptide.getPeptideModSeq());
			}
		}

		HashMap<SearchJobData, TObjectFloatHashMap<String>> retentionTimeMappingsInSeconds=new HashMap<>();
		HashSet<String> addedLibraries=new HashSet<>();
		HashMap<String, CorrelationPeakFrequencyCalculator> ionCounter=new HashMap<String, CorrelationPeakFrequencyCalculator>();

		for (SearchJobData job : jobs) {
			if (job instanceof LibrarySearchJobData) {
				// try reading encyclopedia data directly from results library
				LibrarySearchJobData libjob = (LibrarySearchJobData)job;
				
				File resultLibrary=libjob.getResultLibrary();
				LibraryInterface lib=libjob.getLibrary();
				if (!addedLibraries.contains(lib.getName())) {
					addedLibraries.add(lib.getName());
					Logger.logLine("Adding library ["+lib.getName()+"] to ion tracker...");
					addLibraryToCounter(lib, ionCounter, true, params);
				}
				try {
					LibraryInterface results=BlibToLibraryConverter.getFile(resultLibrary);
					Logger.logLine("Adding results ["+results.getName()+"] to ion tracker...");
					TObjectFloatHashMap<String> rtMapping=addLibraryToCounter(results, ionCounter, false, params);
					retentionTimeMappingsInSeconds.put(job, rtMapping);
					
					
				} catch (EncyclopediaException e) {
					Logger.errorLine("Parsing error indicates "+job.getPercolatorFiles().getPeptideOutputFile().getName()+" isn't from a quantitative search (EncyclopeDIA or XCorDIA):");
					Logger.errorException(e);
				}
			}
		}
		
		int strongAboveThreshold=0;
		HashMap<String,double[]> bestIons=new HashMap<String, double[]>();
		for (Entry<String, CorrelationPeakFrequencyCalculator> entry : ionCounter.entrySet()) {
			String peptideModSeq=entry.getKey();
			
			// always accept if no passing peptides were used
			if (passingPeptideModSeqs.size()==0||passingPeptideModSeqs.contains(peptideModSeq)) {
				double[] ions=entry.getValue().getTopNMasses(params.getEffectiveNumberOfQuantitativePeaks());
				if (ions!=null&&ions.length>0) {
					if (ions.length>=params.getMinNumOfQuantitativePeaks()) {
						strongAboveThreshold++;
					}
					bestIons.put(peptideModSeq, ions);
				}
			}
		}
		Logger.logLine("Found quantitative ions for "+bestIons.size()+" total peptides ("+strongAboveThreshold+" with "+params.getMinNumOfQuantitativePeaks()+" or more high quality peaks) across all runs.");

		return new Pair<HashMap<SearchJobData,TObjectFloatHashMap<String>>, HashMap<String,double[]>>(retentionTimeMappingsInSeconds, bestIons);
	}
	
	private static TObjectFloatHashMap<String> addLibraryToCounter(LibraryInterface lib, HashMap<String, CorrelationPeakFrequencyCalculator> ionCounter, boolean isLibrary, SearchParameters params) {
		TObjectFloatHashMap<String> rtMapping=new TObjectFloatHashMap<String>();
		try {
			ArrayList<LibraryEntry> entries=lib.getAllEntries(false, params.getAAConstants());
			
			for (LibraryEntry entry : entries) {
				String peptideModSeq=entry.getPeptideModSeq();
				rtMapping.put(peptideModSeq, entry.getRetentionTime());
				CorrelationPeakFrequencyCalculator bestIonsMap=ionCounter.get(peptideModSeq);
				if (bestIonsMap==null) {
					bestIonsMap=new CorrelationPeakFrequencyCalculator(params.getFragmentTolerance());
					ionCounter.put(peptideModSeq, bestIonsMap);
				}
				double[] masses=entry.getMassArray();
				float[] intensity=entry.getIntensityArray();
				float[] correlation=entry.getCorrelationArray();
				boolean[] isQuant=entry.getQuantifiedIonsArray();
				for (int i=0; i<correlation.length; i++) {
					boolean passesThreshold = isQuant[i]&&(isLibrary||correlation[i]>=TransitionRefiner.quantitativeCorrelationThreshold);
					float thisCorrelation=0.0f;
					if (correlation[i]>=TransitionRefiner.identificationCorrelationThreshold) {
						thisCorrelation=correlation[i];
					}
					bestIonsMap.increment(masses[i], intensity[i], thisCorrelation, passesThreshold, isLibrary);
				}
			}

		} catch (IOException e) {
			throw new EncyclopediaException("Error parsing results library", e);
		} catch (SQLException e) {
			throw new EncyclopediaException("Error parsing results library", e);
		} catch (DataFormatException e) {
			throw new EncyclopediaException("Error parsing results library", e);
		}
		return rtMapping;
	}
}
