package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.TreeSet;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideQuantExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ChromatogramLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakFrequencyCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class PeakLocationInferrer {
	// alignments are seed (x) to sample (y), in minutes
	private final HashMap<SearchJobData, RetentionTimeFilter> alignmentMap;

	// alignedRTs are as if they were in the seed (x) file
	private final HashMap<String, Float> alignedRTInMinBySequenceMap;

	private final HashMap<String, ChromatogramLibraryEntry> libraryEntryBySequenceMap;
	
	private final HashMap<String, double[]> bestIons;
	private final SearchParameters params;

	PeakLocationInferrer(HashMap<SearchJobData, RetentionTimeFilter> alignmentMap, HashMap<String, ChromatogramLibraryEntry> libraryEntryBySequenceMap, HashMap<String, Float> alignedRTInMinBySequenceMap, HashMap<String, double[]> bestIons, SearchParameters params) {
		this.alignmentMap=alignmentMap;
		this.libraryEntryBySequenceMap=libraryEntryBySequenceMap;
		this.alignedRTInMinBySequenceMap=alignedRTInMinBySequenceMap;
		this.bestIons=bestIons;
		this.params=params;
	}
	
	public Pair<Float, Integer> getTopNIntensity(PeptidePrecursor peptide, TransitionRefinementData data) {
		double[] topN=bestIons.get(peptide.getPeptideModSeq());
		double[] masses=data.getFragmentMassArray();
		float[] correlation=data.getCorrelationArray();
		float[] intensities=data.getIntegrationArray();
		
		if (topN==null||topN.length==0) {
			return data.getTopNIntensity(TransitionRefiner.quantitativeCorrelationThreshold, params.getNumberOfQuantitativePeaks());
		}
		
		float sum=0.0f;
		int added=0;
		for (int i=0; i<topN.length; i++) {
			Optional<Integer> optionalIndex=params.getFragmentTolerance().getIndex(masses, topN[i]);
			if (optionalIndex.isPresent()) {
				int index=optionalIndex.get();
				if (correlation[index]>=TransitionRefiner.translationalQuantitativeCorrelationThreshold) {
					sum+=intensities[index];
					added++;
				}
			}
		}
		return new Pair<Float, Integer>(sum, added);
	}

	public float getRTInSec(SearchJobData job, String peptideModSeq) {
		RetentionTimeFilter f=alignmentMap.get(job);
		Float alignedRTInMin=alignedRTInMinBySequenceMap.get(peptideModSeq);
		if (alignedRTInMin==null) {
			Logger.errorLine("Couldn't find retention time for peptide ("+peptideModSeq+").");
			return -1;
		}
		;
		if (f==null) {
			// job is the seed
			return alignedRTInMin*60f;
		} else {
			return f.getYValue(alignedRTInMin)*60f;
		}
	}

	public static PeakLocationInferrer getAlignmentData(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, ArrayList<PercolatorPeptide> passingPeptides, SearchParameters params) {
		ProgressIndicator subProgress1=new SubProgressIndicator(progress, 0.5f);
		Pair<HashMap<SearchJobData,ArrayList<ChromatogramLibraryEntry>>, HashMap<String,double[]>> pair=getArchetypalPeptides(subProgress1, pecanJobs, passingPeptides, params);
		HashMap<SearchJobData, ArrayList<ChromatogramLibraryEntry>> archetypalPeptides=pair.x;
		HashMap<String, double[]> bestIons=pair.y;

		// get best job
		SearchJobData bestJob=null;
		int max=-1;
		HashMap<String, ChromatogramLibraryEntry> libraryEntryBySequenceMap=new HashMap<String, ChromatogramLibraryEntry>();
		for (Entry<SearchJobData, ArrayList<ChromatogramLibraryEntry>> entry : archetypalPeptides.entrySet()) {
			int length=entry.getValue().size();
			if (length>max) {
				max=length;
				bestJob=entry.getKey();
			}

			// also grab library entry map
			for (ChromatogramLibraryEntry pep : entry.getValue()) {
				libraryEntryBySequenceMap.put(pep.getPeptideModSeq(), pep);
			}
		}

		Logger.logLine("Seed experiment: "+bestJob.getDiaFile().getName());
		Logger.logLine("Seed Percolator file: "+bestJob.getOutputFile().getAbsolutePath());
		ArrayList<PercolatorPeptide> alignmentSeed=PercolatorReader.getPassingPeptidesFromTSV(bestJob.getOutputFile(), bestJob.getParameters().getPercolatorThreshold());
		TObjectFloatHashMap<String> rtsBySequence=new TObjectFloatHashMap<String>();
		for (PercolatorPeptide peptide : alignmentSeed) {
			rtsBySequence.put(peptide.getPeptideModSeq(), peptide.getRT());
		}
		Logger.logLine("Number of anchors in seed file: "+alignmentSeed.size());

		ProgressIndicator subProgress2=new SubProgressIndicator(progress, 0.5f);

		// construct alignments
		HashMap<SearchJobData, RetentionTimeFilter> alignmentMap=new HashMap<SearchJobData, RetentionTimeFilter>();
		HashMap<String, Float> alignedRTInMinBySequenceMap=new HashMap<String, Float>();
		int count=0;
		for (SearchJobData job : pecanJobs) {
			if (job!=bestJob) {
				subProgress2.update(job.getDiaFile().getName()+": RT aligning to seed", count/(float) pecanJobs.size());
				count++;

				ArrayList<XYPoint> points=new ArrayList<XYPoint>();
				ArrayList<PercolatorPeptide> peptides=PercolatorReader.getPassingPeptidesFromTSV(job.getOutputFile(), job.getParameters().getPercolatorThreshold());
				for (PercolatorPeptide peptide : peptides) {
					String seq=peptide.getPeptideModSeq();
					if (rtsBySequence.containsKey(seq)) {
						points.add(new XYPoint(rtsBySequence.get(seq)/60f, peptide.getRT()/60f));
					}
				}
				if (points.size()<10) {
					Logger.errorLine("Not enough points ("+points.size()+" out of align:"+peptides.size()+" and best:"+rtsBySequence.size()+") to compute regression between samples, still trying anyways.");
				}
				
				RetentionTimeFilter alignment=new RetentionTimeFilter(points, bestJob.getDiaFile().getName(), job.getDiaFile().getName());
				alignmentMap.put(job, alignment);
				if (job instanceof EncyclopediaJobData) {
					// try reading encyclopedia data directly from results library
					File resultLibrary=((EncyclopediaJobData) job).getResultLibrary();
					alignment.plot(points, Optional.ofNullable(resultLibrary));
				}

				// align local archetypals to the seed
				ArrayList<ChromatogramLibraryEntry> archetypals=archetypalPeptides.get(job);
				for (ChromatogramLibraryEntry entry : archetypals) {
					float alignedRT=alignment.getXValue(entry.getRetentionTime()/60f);
					alignedRTInMinBySequenceMap.put(entry.getPeptideModSeq(), alignedRT);
				}
			}
		}

		return new PeakLocationInferrer(alignmentMap, libraryEntryBySequenceMap, alignedRTInMinBySequenceMap, bestIons, params);
	}

	/**
	 * divvy up all the globally passing peptides by the individual search that
	 * best identified them
	 * 
	 * @param progress
	 * @param pecanJobs
	 * @param passingPeptides
	 * @return
	 */
	static Pair<HashMap<SearchJobData, ArrayList<ChromatogramLibraryEntry>>, HashMap<String, double[]>> getArchetypalPeptides(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs,
			ArrayList<PercolatorPeptide> passingPeptides, SearchParameters params) {
		int numberOfQuantitativePeaks=params.getNumberOfQuantitativePeaks();
		MassTolerance fragmentTolerance=params.getFragmentTolerance();
		
		// set up data structures
		HashMap<String, SearchJobData> jobsByFile=new HashMap<String, SearchJobData>();
		HashMap<String, ArrayList<PercolatorPeptide>> peptidesByFile=new HashMap<String, ArrayList<PercolatorPeptide>>();
		for (SearchJobData job : pecanJobs) {
			String name=job.getDiaFile().getName();
			name=name.substring(0, name.lastIndexOf('.'));
			jobsByFile.put(name, job);
			peptidesByFile.put(name, new ArrayList<PercolatorPeptide>());
		}

		// the best individual search is imbedded in the psmID
		for (PercolatorPeptide psm : passingPeptides) {
			String name=psm.getFile();
			name=name.substring(0, name.lastIndexOf('.'));
			ArrayList<PercolatorPeptide> list=peptidesByFile.get(name);
			if (list==null) {
				Logger.errorLine("Unexpected file ["+name+"] when parsing Percolator result! Ignoring peptide.");
			} else {
				list.add(psm);
			}
		}

		HashMap<String, PeakFrequencyCalculator> ionCounter=new HashMap<String, PeakFrequencyCalculator>();
		
		// extract out chromatogram library entries
		HashMap<SearchJobData, ArrayList<ChromatogramLibraryEntry>> archetypalPeptides=new HashMap<SearchJobData, ArrayList<ChromatogramLibraryEntry>>();
		float increment=1.0f/pecanJobs.size();
		for (Entry<String, ArrayList<PercolatorPeptide>> entry : peptidesByFile.entrySet()) {
			ProgressIndicator subProgress=new SubProgressIndicator(progress, increment);

			SearchJobData job=jobsByFile.get(entry.getKey());
			
			TreeSet<PeptidePrecursor> peptides=new TreeSet<PeptidePrecursor>(entry.getValue());

			boolean readFromLibraryResult=false;
			if (job instanceof EncyclopediaJobData) {
				// try reading encyclopedia data directly from results library
				File resultLibrary=((EncyclopediaJobData) job).getResultLibrary();
				try {
					LibraryInterface results=BlibToLibraryConverter.getFile(resultLibrary);
					
					/*ArrayList<PeptidePrecursor> recast=new ArrayList<PeptidePrecursor>();
					for (PercolatorPeptide pep : peptides) {
						recast.add(pep);
					}
					Logger.errorLine("Parsed:"+peptides.size());
					HashMap<PeptidePrecursor, ArrayList<LibraryEntry>> entries=results.getEntries(recast, false);
					*/
					
					ArrayList<LibraryEntry> entries=results.getAllEntries(false);
					
					ArrayList<ChromatogramLibraryEntry> bestEntries=new ArrayList<ChromatogramLibraryEntry>();
					for (LibraryEntry libEntry : entries) {
						// all results files are saved as chromatogram libraries
						ChromatogramLibraryEntry chrom=(ChromatogramLibraryEntry)libEntry;
						
						String peptideModSeq=libEntry.getPeptideModSeq();
						PeakFrequencyCalculator bestIonsMap=ionCounter.get(peptideModSeq);
						if (bestIonsMap==null) {
							bestIonsMap=new PeakFrequencyCalculator(fragmentTolerance);
							ionCounter.put(peptideModSeq, bestIonsMap);
						}
						
						double[] masses=chrom.getMassArray();
						float[] correlation=chrom.getCorrelationArray();
						for (int i=0; i<correlation.length; i++) {
							if (correlation[i]>=TransitionRefiner.quantitativeCorrelationThreshold) {
								bestIonsMap.increment(masses[i]);
							}
						}

						if (peptides.contains(libEntry)) {
							bestEntries.add(chrom);
						}
					}
					
					Logger.errorLine(resultLibrary.getName()+"produced Parsed:"+peptides.size()+", BEST:"+bestEntries.size());
					
					archetypalPeptides.put(job, bestEntries);
					
					readFromLibraryResult=true;
				} catch (EncyclopediaException e) {
					Logger.errorLine("Parsing error indicates "+job.getOutputFile().getName()+" isn't from Encyclopedia:");
					Logger.errorException(e);
					readFromLibraryResult=false;
				} catch (IOException e) {
					throw new EncyclopediaException("Error parsing results library", e);
				} catch (SQLException e) {
					throw new EncyclopediaException("Error parsing results library", e);
				} catch (DataFormatException e) {
					throw new EncyclopediaException("Error parsing results library", e);
				}
			}
			
			// if we can't read data from a library result file (e.g. Pecan), then read directly from the DIA file
			if (!readFromLibraryResult) {
				StripeFileInterface stripeFile=MzmlToDIAConverter.getFile(job.getDiaFile(), job.getParameters());
				Logger.logLine("Extracting "+peptides.size()+" Archetypal Peptides from "+job.getDiaFile().getName()+"...");
				subProgress.update(job.getDiaFile().getName()+": Extracting "+peptides.size()+" Archetypal Peptides", 0.00001f);

				LibraryInterface library=null;
				if (job instanceof EncyclopediaJobData) {
					library=((EncyclopediaJobData) job).getLibrary();
				}
				ArrayList<IntegratedLibraryEntry> libraryEntries=PeptideQuantExtractor.parseSearchFeatures(subProgress, job, false, passingPeptides, entry.getValue(), Optional.ofNullable((PeakLocationInferrer)null), stripeFile, library,
						job.getParameters());
				ArrayList<ChromatogramLibraryEntry> recast=new ArrayList<ChromatogramLibraryEntry>();
				for (IntegratedLibraryEntry e : libraryEntries) {
					recast.add(e);
				}
				
				archetypalPeptides.put(job, recast);
				stripeFile.close();
			}
		}
		HashMap<String,double[]> bestIons=new HashMap<String, double[]>();
		for (Entry<String, PeakFrequencyCalculator> entry : ionCounter.entrySet()) {
			String peptideModSeq=entry.getKey();
			double[] ions=entry.getValue().getTopNMasses(numberOfQuantitativePeaks);
			bestIons.put(peptideModSeq, ions);
		}
		return new Pair<HashMap<SearchJobData,ArrayList<ChromatogramLibraryEntry>>, HashMap<String,double[]>>(archetypalPeptides, bestIons);
	}
}
