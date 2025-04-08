package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Optional;
import java.util.StringTokenizer;
import java.util.concurrent.ConcurrentLinkedQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
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
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaLibraryPredictionClient;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaPrecursor;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2019HCDModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt.Prosit2019RTModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
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
		//WQVVDTPGILDHPLEDR	85.748116	2	
		String note="DVVNVLQAVGESLAK	94.66456	3";
		note="SPSDTEGLVK	40.243713	2";
		
		StringTokenizer st=new StringTokenizer(note);
		String peptideModSeq=st.nextToken();
		float retentionTime=Float.parseFloat(st.nextToken())*60f;
		byte precursorCharge=Byte.parseByte(st.nextToken());
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
		System.out.println("Predicted: "+entry.getPeptideModSeq());

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
		System.out.println("Processed: "+poll.getPeptideModSeq());

		float[] expectedIsotopicDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(poll.getPeptideModSeq(), parameters.getAAConstants());
		
		Pair<ArrayList<XYPoint>, ArrayList<XYPoint>[]> precursorTraces=precursors.
				integrateChromatogram(psm.getPrecursorMZ(), psm.getPrecursorCharge(), poll.getRefinementData().getRange(), 
						expectedIsotopicDistribution, parameters.getPrecursorTolerance());
		System.out.println("Found points: "+precursorTraces.getX().size());

		XYTrace precursorTrace = new XYTrace(precursorTraces.x, GraphType.dashedline, "p", Color.red, 1.0f);
		
		//System.out.println(precursorTrace.integrate(data.x.getRange()));
		ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
//		for (int i = 0; i < precursorTraces.y.length; i++) {
//			traces.add(new XYTrace(precursorTraces.y[i], GraphType.boldline, "+"+(i-1)));
//		}
		//System.out.println("Range: "+data.x.getRange()+" vs "+scanStart+" to "+scanStop);
		
		XYTrace fragmentTrace=new XYTrace(poll.getRefinementData().getRtArray().get(), poll.getRefinementData().getMedianChromatogram(), GraphType.dashedline, "f", Color.blue, 1.0f).trim(poll.getRefinementData().getRange());
		System.out.println(fragmentTrace);
		traces.add(precursorTrace.rescaleY(1.0f/(float)precursorTrace.getMaxY()));
		traces.add(fragmentTrace.rescaleY(1.0f/(float)fragmentTrace.getMaxY()));
		

		ArrayList<XYPoint>[] values=XYTrace.alignXYPoints(precursorTrace.getPoints(), fragmentTrace.getPoints());

		XYTrace alignedX = new XYTrace(values[0], GraphType.bighollowpoint, "p2", Color.red, 10.0f);
		XYTrace alignedY = new XYTrace(values[1], GraphType.bighollowpoint, "f2", Color.blue, 10.0f);
        
		traces.add(alignedX.rescaleY(1.0f/(float)alignedX.getMaxY()));
		traces.add(alignedY.rescaleY(1.0f/(float)alignedY.getMaxY()));
		
		System.out.println("Correlation: "+XYTrace.correlate(precursorTrace, fragmentTrace)+" \t"+poll.getRefinementData().getPeptideModSeq()+"\t"+poll.getRefinementData().getApexRT()/60f+"\t"+poll.getRefinementData().getPrecursorCharge());
		
		Charter.launchChart("Retention Time", "Intensity", true, traces.toArray(new XYTrace[0]));
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
