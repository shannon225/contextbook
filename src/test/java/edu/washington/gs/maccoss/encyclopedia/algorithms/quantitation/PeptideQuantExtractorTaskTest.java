package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Optional;
import java.util.concurrent.ConcurrentLinkedQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import junit.framework.TestCase;

public class PeptideQuantExtractorTaskTest extends TestCase {
	public static void main(String[] args) throws Exception {
		File f=new File("/Users/searleb/Documents/encyclopedia/bugs/failed_integration/2020dec03_cobbs_cmv_inf_gpfdia_05_3.dia");
		StripeFile raw=new StripeFile(true);
		raw.openFile(f);
		ConcurrentLinkedQueue<IntegratedLibraryEntry> savedEntries=new ConcurrentLinkedQueue<IntegratedLibraryEntry>();
		
		//2020dec03_cobbs_cmv_inf_gpfdia_05_3.mzML:3876.1084:NLVPMVATVQGQNLK+2
		String peptideModSeq="NLVPMVATVQGQNLK";
		byte precursorCharge=2;
		float retentionTime=3876.1084f;
		int spectrumIndex=103782;
		
		
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		double targetMz=parameters.getAAConstants().getChargedMass(peptideModSeq, precursorCharge);
		float expectedPeakWidth = parameters.getExpectedPeakWidth();
		PSMData psm=new PSMData(new HashSet<>(), spectrumIndex, targetMz, precursorCharge, peptideModSeq, retentionTime, 1, 1, expectedPeakWidth, false, parameters.getAAConstants());
		
		ArrayList<FragmentScan> stripes=raw.getStripes(targetMz, 0, Float.MAX_VALUE, false);
		System.out.println(targetMz+" --> "+stripes.size());
		
		PeptideQuantExtractorTask task=new PeptideQuantExtractorTask(f.getName(), psm, Optional.ofNullable(null), Optional.ofNullable(null), stripes, parameters, savedEntries, false);
		task.process();
		IntegratedLibraryEntry poll=savedEntries.poll();
		System.out.println(poll);
	}
	
	public void testQuantifyPeptide() {
		PSMPeakScorer scorer;
		AnnotatedLibraryEntry unitEntry;
		boolean limitToQuantifiable=false;
		ArrayList<FragmentScan> stripes;
		boolean integrateEverything;
		boolean wasInferred;
		SearchParameters params;
		

		//TransitionRefinementData data = PeptideQuantExtractorTask.quantifyPeptide(scorer, unitEntry, limitToQuantifiable, stripes, integrateEverything, wasInferred, params);
	}
}
