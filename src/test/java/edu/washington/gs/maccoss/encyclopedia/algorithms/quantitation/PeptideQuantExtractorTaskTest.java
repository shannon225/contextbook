package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Optional;
import java.util.concurrent.ConcurrentLinkedQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.LibraryPeakLocationInferrer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.parameters.InstrumentSpecificSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaLibraryPredictionClient;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaPrecursor;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2019HCDModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2020HCDModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2023timsTOFModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.ims.AlphaPeptDeepIMSModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt.ChronologerModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt.Prosit2019RTModel;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;
import junit.framework.TestCase;

public class PeptideQuantExtractorTaskTest extends TestCase {
	public static void main(String[] args) throws Exception {
		final File f=new File("/Users/searleb/Documents/manuscripts/2025/mapms/mapms/2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_01.dia");
		StripeFile raw=new StripeFile(true);
		raw.openFile(f);
		ConcurrentLinkedQueue<IntegratedLibraryEntry> savedEntries=new ConcurrentLinkedQueue<IntegratedLibraryEntry>();
		SearchParameters parameters=InstrumentSpecificSearchParameters.OrbitrapOrbitrap.getDefaultParameters();
		
		//YGIEPTMVVQGVK	74.9156723	2
		String peptideModSeq="YGIEPTMVVQGVK";
		byte precursorCharge=2;
		float retentionTime=74.9156723f*60f;
		int spectrumIndex=0;

		ArrayList<KoinaPrecursor> pred=new ArrayList<KoinaPrecursor>();
		pred.add(new KoinaPrecursor(peptideModSeq, 33f, precursorCharge, parameters.getAAConstants()));
		ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
		models.add(new Prosit2019HCDModel());
		models.add(new Prosit2019RTModel());
		KoinaLibraryPredictionClient client=new KoinaLibraryPredictionClient(models);
		client.generatePredictions(KoinaLibraryPredictionClient.HTTPS_KOINA_WILHELMLAB_ORG_443, pred, new EmptyProgressIndicator(true));
		
		AnnotatedLibraryEntry entry=pred.get(0).toEntry(parameters.getAAConstants(), parameters);
		ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
		entries.add(entry);
		
		HashMap<SearchJobData, TObjectFloatHashMap<String>> rtByPeptideModSeq=new HashMap<SearchJobData, TObjectFloatHashMap<String>>();
		
		LibraryPeakLocationInferrer inferrer=new LibraryPeakLocationInferrer(entries, rtByPeptideModSeq, parameters);
		
		double targetMz=parameters.getAAConstants().getChargedMass(peptideModSeq, precursorCharge);
		float expectedPeakWidth = parameters.getExpectedPeakWidth();
		PSMData psm=new PSMData(new HashSet<>(), spectrumIndex, targetMz, precursorCharge, peptideModSeq, retentionTime, 1, 1, expectedPeakWidth, false, parameters.getAAConstants());
		
		ArrayList<FragmentScan> stripes=raw.getStripes(targetMz, 0, Float.MAX_VALUE, false);
		System.out.println("mz:"+targetMz+", rt: "+38.59886932f*60f);

		Logger.logLine("Processing precursors scans...");
		PrecursorScanMap precursors=new PrecursorScanMap(raw.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));
		
		PeptideQuantExtractorTask task=new PeptideQuantExtractorTask(f.getName(), psm, Optional.ofNullable(inferrer), Optional.ofNullable(null), stripes, Optional.ofNullable(precursors), parameters, savedEntries, false);
		task.process();
		IntegratedLibraryEntry poll=savedEntries.poll();
		System.out.println(poll.getPeptideModSeq());
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
