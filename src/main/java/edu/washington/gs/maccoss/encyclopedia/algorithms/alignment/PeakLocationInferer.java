package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideQuantExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class PeakLocationInferer {
	// alignments are seed (x) to sample (y)
	HashMap<SearchJobData, Function> alignmentMap;
	
	// alignedRTs are as if they were in the seed (x) file
	HashMap<String, Float> alignedRTInSecBySequenceMap;
	
	HashMap<String, IntegratedLibraryEntry> libraryEntryBySequenceMap;
	
	PeakLocationInferer(HashMap<SearchJobData, Function> alignmentMap, HashMap<String, IntegratedLibraryEntry> libraryEntryBySequenceMap, HashMap<String, Float> alignedRTInSecBySequenceMap) {
		this.alignmentMap=alignmentMap;
		this.libraryEntryBySequenceMap=libraryEntryBySequenceMap;
		this.alignedRTInSecBySequenceMap=alignedRTInSecBySequenceMap;
	}
	
	public float getRTInSec(SearchJobData job, String peptide) {
		Function f=alignmentMap.get(job);
		Float alignedRT=alignedRTInSecBySequenceMap.get(peptide);
		if (alignedRT==null) {
			Logger.errorLine("Couldn't find retention time for peptide ("+peptide+").");
			return -1;
		};
		if (f==null) {
			// job is the seed 
			return alignedRT;
		} else {
			return f.getYValue(alignedRT);
		}
	}

	public static PeakLocationInferer getAlignmentData(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, File blibFile, ArrayList<PercolatorPeptide> passingPeptides) {
		ProgressIndicator subProgress1=new SubProgressIndicator(progress, 0.5f);
		HashMap<SearchJobData, ArrayList<IntegratedLibraryEntry>> archetypalPeptides=getArchetypalPeptides(subProgress1, pecanJobs, blibFile, passingPeptides);
		
		// get best job
		SearchJobData bestJob=null;
		int max=-1;
		HashMap<String, IntegratedLibraryEntry> libraryEntryBySequenceMap=new HashMap<String, IntegratedLibraryEntry>();
		for (Entry<SearchJobData, ArrayList<IntegratedLibraryEntry>> entry : archetypalPeptides.entrySet()) {
			int length=entry.getValue().size();
			if (length>max) {
				max=length;
				bestJob=entry.getKey();
			}
			
			// also grab library entry map
			for (IntegratedLibraryEntry pep : entry.getValue()) {
				libraryEntryBySequenceMap.put(pep.getPeptideModSeq(), pep);
			}
		}
		
		ArrayList<PercolatorPeptide> alignmentSeed=PercolatorReader.getPassingPeptidesFromTSV(bestJob.getOutputFile(), bestJob.getParameters().getPercolatorThreshold());
		TObjectFloatHashMap<String> rtsBySequence=new TObjectFloatHashMap<String>();
		for (PercolatorPeptide peptide : alignmentSeed) {
			rtsBySequence.put(peptide.getPeptideSequence(), peptide.getRT());
		}

		ProgressIndicator subProgress2=new SubProgressIndicator(progress, 0.5f);

		// construct alignments
		HashMap<SearchJobData, Function> alignmentMap=new HashMap<SearchJobData, Function>();
		HashMap<String, Float> alignedRTInSecBySequenceMap=new HashMap<String, Float>();
		int count=0;
		for (SearchJobData job : pecanJobs) {
			if (job!=bestJob) {
				subProgress2.update(job.getDiaFile().getName()+": RT aligning to seed", count/(float)pecanJobs.size());
				count++;
				
				ArrayList<XYPoint> points=new ArrayList<XYPoint>();
				ArrayList<PercolatorPeptide> peptides=PercolatorReader.getPassingPeptidesFromTSV(job.getOutputFile(), job.getParameters().getPercolatorThreshold());
				for (PercolatorPeptide peptide : peptides) {
					String seq=peptide.getPeptideSequence();
					if (rtsBySequence.containsKey(seq)) {
						points.add(new XYPoint(rtsBySequence.get(seq), peptide.getRT()));
					}
				}
				if (points.size()<10) {
					Logger.errorLine("Not enough points ("+points.size()+") to compute regression between samples, still trying anyways.");
				}
				TwoDimensionalKDE twoDimKDE=new TwoDimensionalKDE(points);
				Function alignment=twoDimKDE.trace();
				alignmentMap.put(job, alignment);
				
				// align local archetyals to the seed
				ArrayList<IntegratedLibraryEntry> archetypals=archetypalPeptides.get(job);
				for (IntegratedLibraryEntry entry : archetypals) {
					float alignedRT=alignment.getXValue(entry.getRetentionTime());
					alignedRTInSecBySequenceMap.put(entry.getPeptideModSeq(), alignedRT);
				}
			}
		}
		
		return new PeakLocationInferer(alignmentMap, libraryEntryBySequenceMap, alignedRTInSecBySequenceMap);
	}
	
	/**
	 * divvy up all the globally passing peptides by the individual search that best identified them 
	 * @param progress
	 * @param pecanJobs
	 * @param blibFile
	 * @param passingPeptides
	 * @return
	 */
	static HashMap<SearchJobData, ArrayList<IntegratedLibraryEntry>> getArchetypalPeptides(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, File blibFile, ArrayList<PercolatorPeptide> passingPeptides) {
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
		HashMap<SearchJobData, ArrayList<IntegratedLibraryEntry>> archetypalPeptides=new HashMap<SearchJobData, ArrayList<IntegratedLibraryEntry>>();
		float increment=1.0f/pecanJobs.size();
		for (Entry<String, ArrayList<PercolatorPeptide>> entry : peptidesByFile.entrySet()) {
			ProgressIndicator subProgress=new SubProgressIndicator(progress, increment);
			
			SearchJobData job=jobsByFile.get(entry.getKey());
			ArrayList<PercolatorPeptide> peptides=entry.getValue();
			
			StripeFileInterface stripeFile=MzmlToDIAConverter.getFile(job.getDiaFile(), job.getParameters());
			Logger.logLine("Extracting "+peptides.size()+" Archetypal Peptides from "+job.getDiaFile().getName()+"...");
			subProgress.update(job.getDiaFile().getName()+": Extracting "+peptides.size()+" Archetypal Peptides", 0.00001f);

			LibraryInterface library=null;
			if (job instanceof EncyclopediaJobData) {
				library=((EncyclopediaJobData)job).getLibrary();
			}
			ArrayList<IntegratedLibraryEntry> libraryEntries=PeptideQuantExtractor.parseSearchFeatures(subProgress, job.getFeatureFile(), false, passingPeptides, peptides, stripeFile, library, job.getParameters());
			archetypalPeptides.put(job, libraryEntries);
			stripeFile.close();
		}
		return archetypalPeptides;
	}
}
