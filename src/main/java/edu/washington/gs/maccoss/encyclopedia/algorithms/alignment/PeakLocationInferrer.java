package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Optional;
import java.util.Map.Entry;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideQuantExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ChromatogramLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class PeakLocationInferrer {
	// alignments are seed (x) to sample (y), in minutes
	HashMap<SearchJobData, RetentionTimeFilter> alignmentMap;

	// alignedRTs are as if they were in the seed (x) file
	HashMap<String, Float> alignedRTInMinBySequenceMap;

	HashMap<String, ChromatogramLibraryEntry> libraryEntryBySequenceMap;

	PeakLocationInferrer(HashMap<SearchJobData, RetentionTimeFilter> alignmentMap, HashMap<String, ChromatogramLibraryEntry> libraryEntryBySequenceMap, HashMap<String, Float> alignedRTInMinBySequenceMap) {
		this.alignmentMap=alignmentMap;
		this.libraryEntryBySequenceMap=libraryEntryBySequenceMap;
		this.alignedRTInMinBySequenceMap=alignedRTInMinBySequenceMap;
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

	public static PeakLocationInferrer getAlignmentData(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, ArrayList<PercolatorPeptide> passingPeptides) {
		ProgressIndicator subProgress1=new SubProgressIndicator(progress, 0.5f);
		HashMap<SearchJobData, ArrayList<ChromatogramLibraryEntry>> archetypalPeptides=getArchetypalPeptides(subProgress1, pecanJobs, passingPeptides);

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

		ArrayList<PercolatorPeptide> alignmentSeed=PercolatorReader.getPassingPeptidesFromTSV(bestJob.getOutputFile(), bestJob.getParameters().getPercolatorThreshold());
		TObjectFloatHashMap<String> rtsBySequence=new TObjectFloatHashMap<String>();
		for (PercolatorPeptide peptide : alignmentSeed) {
			rtsBySequence.put(peptide.getPeptideModSeq(), peptide.getRT());
		}

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

				// align local archetyals to the seed
				ArrayList<ChromatogramLibraryEntry> archetypals=archetypalPeptides.get(job);
				for (ChromatogramLibraryEntry entry : archetypals) {
					float alignedRT=alignment.getXValue(entry.getRetentionTime()/60f);
					alignedRTInMinBySequenceMap.put(entry.getPeptideModSeq(), alignedRT);
				}
			}
		}

		return new PeakLocationInferrer(alignmentMap, libraryEntryBySequenceMap, alignedRTInMinBySequenceMap);
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
	static HashMap<SearchJobData, ArrayList<ChromatogramLibraryEntry>> getArchetypalPeptides(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs,
			ArrayList<PercolatorPeptide> passingPeptides) {
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

		// extract out chromatogram library entries
		HashMap<SearchJobData, ArrayList<ChromatogramLibraryEntry>> archetypalPeptides=new HashMap<SearchJobData, ArrayList<ChromatogramLibraryEntry>>();
		float increment=1.0f/pecanJobs.size();
		for (Entry<String, ArrayList<PercolatorPeptide>> entry : peptidesByFile.entrySet()) {
			ProgressIndicator subProgress=new SubProgressIndicator(progress, increment);

			SearchJobData job=jobsByFile.get(entry.getKey());
			ArrayList<PercolatorPeptide> peptides=entry.getValue();

			boolean readFromLibraryResult=false;
			if (job instanceof EncyclopediaJobData) {
				// try reading encyclopedia data directly from results library
				File resultLibrary=((EncyclopediaJobData) job).getResultLibrary();
				try {
					LibraryInterface results=BlibToLibraryConverter.getFile(resultLibrary);
					
					ArrayList<PeptidePrecursor> recast=new ArrayList<PeptidePrecursor>();
					for (PercolatorPeptide pep : peptides) {
						recast.add(pep);
					}
					System.out.println("Parsed:"+peptides.size());
					HashMap<PeptidePrecursor, ArrayList<LibraryEntry>> entries=results.getEntries(recast, false);
					
					// all results files are saved as chromatogram libraries
					ArrayList<ChromatogramLibraryEntry> bestEntries=new ArrayList<ChromatogramLibraryEntry>();
					for (ArrayList<LibraryEntry> resultEntries : entries.values()) {
						if (resultEntries.size()>0) {
							bestEntries.add((ChromatogramLibraryEntry)resultEntries.get(0));
						}
					}
					System.out.println("BEST:"+bestEntries.size());
					
					archetypalPeptides.put(job, bestEntries);
					
					readFromLibraryResult=true;
				} catch (EncyclopediaException e) {
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
				ArrayList<IntegratedLibraryEntry> libraryEntries=PeptideQuantExtractor.parseSearchFeatures(subProgress, job, false, passingPeptides, peptides, Optional.ofNullable((PeakLocationInferrer)null), stripeFile, library,
						job.getParameters());
				ArrayList<ChromatogramLibraryEntry> recast=new ArrayList<ChromatogramLibraryEntry>();
				for (IntegratedLibraryEntry e : libraryEntries) {
					recast.add(e);
				}
				
				archetypalPeptides.put(job, recast);
				stripeFile.close();
			}
		}
		return archetypalPeptides;
	}
}
