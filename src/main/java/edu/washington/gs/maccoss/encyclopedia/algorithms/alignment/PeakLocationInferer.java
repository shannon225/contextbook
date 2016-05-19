package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideQuantExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;

public class PeakLocationInferer {
	HashMap<SearchJobData, RetentionTimeFilter> alignmentMap;
	HashMap<IntegratedLibraryEntry, Float> alignedRTInSecMap;
	
	/*static PeakLocationInferer getAlignmentData(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, File blibFile, ArrayList<ScoredObject<String>> passingPeptides) {
		ProgressIndicator subProgress1=new SubProgressIndicator(progress, 0.5f);
		HashMap<SearchJobData, ArrayList<IntegratedLibraryEntry>> archetypalPeptides=getArchetypalPeptides(subProgress1, pecanJobs, blibFile, passingPeptides);
		
		SearchJobData bestJob=null;
		int max=-1;
		for (Entry<SearchJobData, ArrayList<IntegratedLibraryEntry>> entry : archetypalPeptides.entrySet()) {
			int length=entry.getValue().size();
			if (length>max) {
				max=length;
				bestJob=entry.getKey();
			}
		}
		
		ArrayList<ScoredObject<String>> 

		ProgressIndicator subProgress2=new SubProgressIndicator(progress, 0.5f);

		for (SearchJobData job : pecanJobs) {
			ArrayList<ScoredObject<String>> peptides=PercolatorReader.getPassingPeptidesFromTSV(job.getOutputFile(), job.getParameters().getPercolatorThreshold());
			for (ScoredObject<String> peptide : peptides) {
				
			}
		}
	}*/
	
	static HashMap<SearchJobData, ArrayList<IntegratedLibraryEntry>> getArchetypalPeptides(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, File blibFile, ArrayList<PercolatorPeptide> passingPeptides) {
		HashMap<String, SearchJobData> jobsByFile=new HashMap<String, SearchJobData>();
		HashMap<String, ArrayList<PercolatorPeptide>> peptidesByFile=new HashMap<String, ArrayList<PercolatorPeptide>>();
		for (SearchJobData job : pecanJobs) {
			String name=job.getDiaFile().getName();
			name=name.substring(0, name.lastIndexOf('.'));
			jobsByFile.put(name, job);
			peptidesByFile.put(name, new ArrayList<PercolatorPeptide>());
		}

		for (PercolatorPeptide psm : passingPeptides) {
			String name=PercolatorReader.getFile(psm.getPsmID());
			name=name.substring(0, name.lastIndexOf('.'));
			ArrayList<PercolatorPeptide> list=peptidesByFile.get(name);
			if (list==null) {
				Logger.errorLine("Unexpected file ["+name+"] when parsing Percolator result! Ignoring peptide.");
			} else {
				list.add(psm);
			}
		}

		HashMap<SearchJobData, ArrayList<IntegratedLibraryEntry>> archetypalPeptides=new HashMap<SearchJobData, ArrayList<IntegratedLibraryEntry>>();
		float increment=1.0f/pecanJobs.size();
		for (Entry<String, ArrayList<PercolatorPeptide>> entry : peptidesByFile.entrySet()) {
			ProgressIndicator subProgress=new SubProgressIndicator(progress, increment);
			
			SearchJobData job=jobsByFile.get(entry.getKey());
			ArrayList<PercolatorPeptide> peptides=entry.getValue();
			
			StripeFileInterface stripeFile=MzmlToDIAConverter.getFile(job.getDiaFile(), job.getParameters());
			Logger.logLine("Extracting "+peptides.size()+" Archetypal Peptides from "+job.getDiaFile().getName()+"...");
			subProgress.update(job.getDiaFile().getName()+": Extracting "+peptides.size()+" Archetypal Peptides", 0.00001f);

			ArrayList<IntegratedLibraryEntry> libraryEntries=PeptideQuantExtractor.parseSearchFeatures(subProgress, job.getFeatureFile(), false, passingPeptides, peptides, stripeFile, Optional.ofNullable((LibraryInterface)null), job.getParameters());
			archetypalPeptides.put(job, libraryEntries);
			stripeFile.close();
		}
		return archetypalPeptides;
	}
}
