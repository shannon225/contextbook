package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.TreeMap;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ChromatogramLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class AlternatePeakLocationInferrer {
	public static PeakLocationInferrerInterface getAlignmentData(ProgressIndicator progress, List<? extends SearchJobData> pecanJobs, ArrayList<PercolatorPeptide> passingPeptides, SearchParameters params) {
		ProgressIndicator subProgress1=new SubProgressIndicator(progress, 0.5f);
		Pair<HashMap<SearchJobData,TObjectFloatHashMap<String>>, HashMap<String,double[]>> pair=getArchetypals(subProgress1, pecanJobs, passingPeptides, params);
		return getInferrer(progress, pecanJobs, pair, params);
	}

	public static PeakLocationInferrerInterface getInferrer(ProgressIndicator progress,
			List<? extends SearchJobData> pecanJobs,
			Pair<HashMap<SearchJobData, TObjectFloatHashMap<String>>, HashMap<String, double[]>> pair,
			SearchParameters params) {
		HashMap<SearchJobData, TObjectFloatHashMap<String>> peptideMappings=pair.x;
		HashMap<String, double[]> bestIons=pair.y;

		// get best job
		SearchJobData bestJob=null;
		int max=-1;
		for (Entry<SearchJobData, TObjectFloatHashMap<String>> entry : peptideMappings.entrySet()) {
			int length=entry.getValue().size();
			if (length>max) {
				max=length;
				bestJob=entry.getKey();
			} else if (length==max&&bestJob!=null) {
				if (bestJob.getDiaFileReader().getOriginalFileName().compareTo(entry.getKey().getDiaFileReader().getOriginalFileName())>0) {
					bestJob=entry.getKey();
				}
			}
		}
		Logger.logLine("Setting "+bestJob.getDiaFileReader().getOriginalFileName()+" as the seed experiment.");
		TObjectFloatHashMap<String> bestRTInSec=peptideMappings.get(bestJob);

		// construct alignments
		HashMap<SearchJobData, RetentionTimeAlignmentInterface> alignmentMap=new HashMap<SearchJobData, RetentionTimeAlignmentInterface>();
		HashMap<SearchJobData, List<RetentionTimeAlignmentInterface.AlignmentDataPoint>> alignmentDataMap = new HashMap<>();
		HashMap<String, Float> alignedRTInMinBySequenceMap=new HashMap<String, Float>();
		// add all bestJob archetypals
		TObjectFloatHashMap<String> archetypals=peptideMappings.get(bestJob);
		archetypals.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String a, float b) {
				alignedRTInMinBySequenceMap.put(a, b/60.0f);
				return true;
			}
		});
		
		ProgressIndicator subProgress2=new SubProgressIndicator(progress, 0.5f);
		int count=0;
		for (SearchJobData job : pecanJobs) {
			if (job!=bestJob) {
				subProgress2.update(job.getDiaFileReader().getOriginalFileName()+": RT aligning to seed", count/(float) pecanJobs.size());
				count++;

				TObjectFloatHashMap<String> rtInSec=peptideMappings.get(job);
				ArrayList<XYPoint> points=new ArrayList<XYPoint>();
				
				bestRTInSec.forEachEntry(new TObjectFloatProcedure<String>() {
					@Override
					public boolean execute(String a, float b) {
						float alt=rtInSec.get(a);
						if (rtInSec.getNoEntryValue()!=alt) {
							points.add(new XYPoint(b/60f, alt/60f)); // both in minutes
						}
						return true;
					}
				});
				
				if (points.size()<10) {
					Logger.errorLine("Not enough points ("+points.size()+" out of align:"+rtInSec.size()+" and best:"+bestRTInSec.size()+") to compute regression between samples, still trying anyways.");
				}
				
				RetentionTimeAlignmentInterface alignment=RetentionTimeFilter.getFilter(points, bestJob.getDiaFileReader().getOriginalFileName(), job.getDiaFileReader().getOriginalFileName());
				alignmentMap.put(job, alignment);

				File saveFileSeed = new File(job.getPercolatorFiles().getPeptideOutputFile().getParentFile(), job.getDiaFileReader().getOriginalFileName());

				final List<RetentionTimeAlignmentInterface.AlignmentDataPoint> alignmentResults =
						alignment.plot(points, Optional.ofNullable(saveFileSeed));
				alignmentDataMap.put(job, alignmentResults);

				// align local archetypals to the seed
				rtInSec.forEachEntry(new TObjectFloatProcedure<String>() {
					@Override
					public boolean execute(String a, float b) {
						float alignedRT=alignment.getXValue(b/60f);
						alignedRTInMinBySequenceMap.put(a, alignedRT);
						return true;
					}
				});
			}
		}

		return new SimplePeakLocationInferrer(alignmentMap, alignmentDataMap, alignedRTInMinBySequenceMap, bestIons, params);
	}

	static Pair<HashMap<SearchJobData, TObjectFloatHashMap<String>>, HashMap<String, double[]>> getArchetypals(ProgressIndicator progress, List<? extends SearchJobData> jobs, ArrayList<PercolatorPeptide> passingPeptides, SearchParameters params) {
		int numberOfQuantitativePeaks=params.getEffectiveNumberOfQuantitativePeaks();
		MassTolerance fragmentTolerance=params.getFragmentTolerance();

		HashMap<SearchJobData, TObjectFloatHashMap<String>> retentionTimeMappingsInSeconds=new HashMap<>();
		HashMap<String, CorrelationPeakFrequencyCalculator> ionCounter=new HashMap<String, CorrelationPeakFrequencyCalculator>();
		HashMap<String, CorrelationPeakFrequencyCalculator> weakIonCounter=new HashMap<String, CorrelationPeakFrequencyCalculator>();
		
		for (SearchJobData job : jobs) {
			if (job instanceof QuantitativeSearchJobData) {
				TObjectFloatHashMap<String> rtMapping=new TObjectFloatHashMap<>();
				retentionTimeMappingsInSeconds.put(job, rtMapping);
				// try reading encyclopedia data directly from results library
				File resultLibrary=((QuantitativeSearchJobData) job).getResultLibrary();
				try {
					LibraryInterface results=BlibToLibraryConverter.getFile(resultLibrary);
					ArrayList<LibraryEntry> entries=results.getAllEntries(false, params.getAAConstants());

					TreeMap<PeptidePrecursor, LibraryEntry> fastLookupPeptides=new TreeMap<PeptidePrecursor, LibraryEntry>();
					for (LibraryEntry libraryEntry : entries) {
						fastLookupPeptides.put(libraryEntry, libraryEntry);
					}
					
					ArrayList<PercolatorPeptide> missingPeptides=new ArrayList<PercolatorPeptide>();
					for (PercolatorPeptide peptide : passingPeptides) {
						LibraryEntry libEntry=fastLookupPeptides.get(peptide);
						if (libEntry==null||!(libEntry instanceof ChromatogramLibraryEntry)) {
							missingPeptides.add(peptide);
							continue;
						}
						// all results files are saved as chromatogram libraries
						ChromatogramLibraryEntry chrom=(ChromatogramLibraryEntry)libEntry;
						rtMapping.put(peptide.getPeptideModSeq(), chrom.getRetentionTime());
						
						String peptideModSeq=libEntry.getPeptideModSeq();
						CorrelationPeakFrequencyCalculator bestIonsMap=ionCounter.get(peptideModSeq);
						if (bestIonsMap==null) {
							bestIonsMap=new CorrelationPeakFrequencyCalculator(fragmentTolerance);
							ionCounter.put(peptideModSeq, bestIonsMap);
						}
						CorrelationPeakFrequencyCalculator bestWeakIonsMap=weakIonCounter.get(peptideModSeq);
						if (bestWeakIonsMap==null) { // note, weak ions aren't used anymore
							bestWeakIonsMap=new CorrelationPeakFrequencyCalculator(fragmentTolerance);
							weakIonCounter.put(peptideModSeq, bestWeakIonsMap);
						}
						double[] masses=chrom.getMassArray();
						float[] intensity=chrom.getIntensityArray();
						float[] correlation=chrom.getCorrelationArray();
						for (int i=0; i<correlation.length; i++) {
							bestIonsMap.increment(masses[i], intensity[i], correlation[i], correlation[i]>=TransitionRefiner.quantitativeCorrelationThreshold);
							bestWeakIonsMap.increment(masses[i], intensity[i], correlation[i], correlation[i]>=TransitionRefiner.identificationCorrelationThreshold);
						}
					}
					
					if (missingPeptides.size()>0) {
						int modified=0;
						for (PercolatorPeptide p : missingPeptides) {
							if (p.getPeptideModSeq().indexOf('[')>=0) {
								modified++;
							}
						}
						Logger.logLine("Couldn't collect global quantitative ion data from this specific sample on "+missingPeptides.size()+" of "+passingPeptides.size()+", (where "+modified+" were modified) at the specified FDR threshold");
					}
					
				} catch (EncyclopediaException e) {
					Logger.errorLine("Parsing error indicates "+job.getPercolatorFiles().getPeptideOutputFile().getName()+" isn't from a quantitative search (EncyclopeDIA or XCorDIA):");
					Logger.errorException(e);
				} catch (IOException e) {
					throw new EncyclopediaException("Error parsing results library", e);
				} catch (SQLException e) {
					throw new EncyclopediaException("Error parsing results library", e);
				} catch (DataFormatException e) {
					throw new EncyclopediaException("Error parsing results library", e);
				}
			}
		}
		
		int strongAboveThreshold=0;
		HashMap<String,double[]> bestIons=new HashMap<String, double[]>();
		for (Entry<String, CorrelationPeakFrequencyCalculator> entry : ionCounter.entrySet()) {
			String peptideModSeq=entry.getKey();
			double[] ions=entry.getValue().getTopNMasses(numberOfQuantitativePeaks);
			if (ions==null||ions.length==0) {
				//double[] altIons=weakIonCounter.get(peptideModSeq).getTopNMasses(numberOfQuantitativePeaks);
				//bestIons.put(peptideModSeq, altIons);
			} else {
				if (ions.length>=params.getMinNumOfQuantitativePeaks()) {
					strongAboveThreshold++;
				}
				bestIons.put(peptideModSeq, ions);
			}
		}
		Logger.logLine("Found quantitative ions for "+bestIons.size()+" total peptides ("+strongAboveThreshold+" with "+params.getMinNumOfQuantitativePeaks()+" or more high quality peaks) across all runs.");

		return new Pair<HashMap<SearchJobData,TObjectFloatHashMap<String>>, HashMap<String,double[]>>(retentionTimeMappingsInSeconds, bestIons);
	}
}
