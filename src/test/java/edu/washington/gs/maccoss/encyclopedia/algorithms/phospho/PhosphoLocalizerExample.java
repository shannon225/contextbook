package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.awt.Color;
import java.awt.Dimension;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class PhosphoLocalizerExample {

	public static void main(String[] args) throws Exception {
		File libraryFile=new File("/Users/searleb/Documents/school/localization_manuscript/hela_phospho/VillenJ_Exactive_HumanPhosphoproteome.elib");
		File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/phospho_repeats/20170430_HeLa_phosp_DIA_B_01_170506220515.dia");

		LibraryFile library=new LibraryFile();
		library.openFile(libraryFile);
		
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		UnitBackgroundFrequencyCalculator background=new UnitBackgroundFrequencyCalculator(0.01f);
		float duration=stripefile.getGradientLength()/20.0f;
		
		PhosphoLocalizer localizer=new PhosphoLocalizer(stripefile, background, parameters);
		
		String peptideModSeq="RPMEEDGEEKS[+80.0]PSK";
		float retentionTime=1434.3873f;
		byte precursorCharge=3;
		/*
		String peptideModSeq="AVT[+80.0]PVPTKTEEVSNLK";
		float retentionTime=3658.4482f;
		byte precursorCharge=3;
		*/
		/*
		String peptideModSeq="IDDRDS[+80.0]DEEGASDR";
		float retentionTime=1606f;
		byte precursorCharge=2;
		*/
		/*
		String peptideModSeq="RAGDLLEDS[+80.0]PKRPK";
		float retentionTime=36*60f;
		byte precursorCharge=3;
		*/
		
		LibraryEntry libentry=library.getEntries(peptideModSeq, precursorCharge, false).get(0);
		double precursorMz=parameters.getAAConstants().getChargedMass(peptideModSeq, precursorCharge);
		
		ArrayList<Stripe> stripes=stripefile.getStripes(precursorMz, 0, Float.MAX_VALUE, false);
		ArrayList<String> permutations=PhosphoPermuter.getPermutations(peptideModSeq, parameters.getAAConstants());
		PhosphoLocalizationData phosphoData=localizer.extractPhosphoFormsFromStripes(peptideModSeq, precursorMz, precursorCharge, permutations, retentionTime, stripes, true);

		System.out.println("Just off of localization ions");
		ArrayList<String> keys=new ArrayList<String>(phosphoData.getPassingForms().keySet());
		for (String sequenceKey : keys) {

			XYPoint point=phosphoData.getLocalizationScores().get(sequenceKey);
			float rt=(float)point.x;
			float localizationScore=(float)point.y;
			
			System.out.println(sequenceKey+"\t"+phosphoData.getPassingForms().get(sequenceKey).getApexRT()+"\t"+rt+"\t"+localizationScore);
		}

		EncyclopediaOneScorer scorer=new EncyclopediaOneScorer(parameters, background);
		FragmentationModel model=new FragmentationModel(libentry.getPeptideModSeq(), parameters.getAAConstants());
		FragmentIon[] ions=model.getPrimaryIonObjects(parameters.getFragType(), libentry.getPrecursorCharge());
		TFloatFloatHashMap primary=new TFloatFloatHashMap();
		for (int i=0; i<stripes.size(); i++) {
			Stripe stripe=stripes.get(i);
			if (stripe.getScanStartTime()>retentionTime-duration&&stripe.getScanStartTime()<retentionTime+duration) {
				primary.put(stripe.getScanStartTime()/60f, scorer.score(libentry, stripe, ions));
			}
		}
		
		HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> allVsUniqueList=phosphoData.getScoreTraces();
		ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
		for (Entry<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> entry : allVsUniqueList.entrySet()) {
			String seq=entry.getKey();
			Pair<TFloatFloatHashMap, TFloatFloatHashMap> pair=entry.getValue();
			Color color=RandomGenerator.randomColor(seq.hashCode());
			//traces.add(new XYTrace(pair.x, GraphType.line, "ALL_"+seq, color, 5.0f));
			traces.add(new XYTrace(pair.y, GraphType.line, "UNI_"+seq, color, 3.0f));
		}
		traces.add(new XYTrace(primary, GraphType.boldline, "primary"));
		
		Charter.launchChart("Retention Time (All Ions)", "Score", true, new Dimension(1000, 400), traces.toArray(new XYTrace[traces.size()]));

		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));
		Range range=null;
		for (Range thisRange : stripefile.getRanges().keySet()) {
			if (thisRange.contains(precursorMz)) {
				range=thisRange;
				break;
			}
		}
		float dutyCycle=stripefile.getRanges().get(range);
		ArrayList<LibraryEntry> entries=new ArrayList<>();
		entries.add(libentry);
		BlockingQueue<PeptideScoringResult> resultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
		CAPSiLOneScoringTask task=new CAPSiLOneScoringTask(scorer, entries, stripes, dutyCycle, precursors, localizer, CAPSiLScoringBreadthType.RECALIBRATED_PEAK_WIDTH, resultsQueue, parameters);
		task.call();

		System.out.println("Based on all ions");
		int index=0;
		while (!resultsQueue.isEmpty()) {
			if (!resultsQueue.isEmpty()) {
				PeptideScoringResult result=resultsQueue.take();
				ArrayList<Pair<ScoredObject<Stripe>, float[]>> data=result.getGoodStripes();
				index++;
				for (Pair<ScoredObject<Stripe>, float[]> pair : data) {
					System.out.println(index+") "+result.getEntry().getPeptideModSeq()+"\t"+pair.x.x+"\t"+pair.x.y.getScanStartTime());
				}
			} else {
				Thread.sleep(10);
			}
		}
	}
}
