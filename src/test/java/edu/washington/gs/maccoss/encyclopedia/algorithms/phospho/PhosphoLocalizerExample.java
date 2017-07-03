package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.awt.Color;
import java.awt.Dimension;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;
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
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class PhosphoLocalizerExample {

	@SuppressWarnings("unused")
	public static void main(String[] args) throws Exception {
		File libraryFile=new File("/Users/searleb/Documents/school/localization_manuscript/hela_phospho/VillenJ_Exactive_HumanPhosphoproteome.elib");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/22jun2016_mcf7_phospho_1b.dia");
		File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/22jun2016_mcf7_phospho_1a.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/hela_phospho/110515_bcs_hela_phospho_starved_20mz_500_900.dia");

		//File libraryFile=new File("/Users/searleb/Documents/phospho_localization/data/VillenJ_Exactive_HumanPhosphoproteome.elib");
		//File diaFile=new File("/Volumes/BriansSSD/phospho/hela_repeats_prism/thesaurus_recalibrated_20p/20170430_HeLa_phosp_DIA_B_04.dia");
		
		//File libraryFile=new File("/Users/searleb/Documents/phospho_localization/data/VillenJ_Exactive_HumanPhosphoproteome.elib");
		//File diaFile=new File("/Users/searleb/Documents/phospho_localization/data/hela/110515_bcs_hela_phospho_starved_20mz_500_900.dia");

		LibraryFile library=new LibraryFile();
		library.openFile(libraryFile);
		
		HashMap<String, String> defaults=SearchParameterParser.getDefaultParameters();
		defaults.put("-localizationModification", "Phosphorylation");
		defaults.put("-ptol", "16.67");
		defaults.put("-ftol", "16.67");
		defaults.put("-lftol", "16.67");
		defaults.put("-scoringBreadthType", "uncal20");
		
		SearchParameters parameters=SearchParameterParser.parseParameters(defaults);
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		UnitBackgroundFrequencyCalculator background=new UnitBackgroundFrequencyCalculator(0.01f);
		float duration=stripefile.getGradientLength()/20.0f;
		
		PhosphoLocalizer localizer=new PhosphoLocalizer(stripefile, PeptideModification.phosphorylation, background, parameters);
		
		String peptideModSeq;
		float retentionTime;
		byte precursorCharge;
		if (false) {
			// from 110515_bcs_hela_phospho_starved_20mz_500_900
			peptideModSeq="RPMEEDGEEKS[+80.0]PSK";
			retentionTime=32*60f;
			precursorCharge=3;
		} else if (false) {
			// repeat 4
			peptideModSeq="ATAPQTQHVSPMR";
			retentionTime=2250.7158203125f;
			precursorCharge=3;
		} else if (false) {
			// repeat 4
			peptideModSeq="LGIAVIHGEAQDAESDLVDGRHS[+80.0]PPMVR";
			retentionTime=5055.13623046875f;
			precursorCharge=4;
		} else if (false) {
			// repeat 4
			peptideModSeq="GRPPAEKLS[+80.0]PNPPNLTK";
			retentionTime=3180.820922851562f;
			precursorCharge=3;
		} else if (false) {
			peptideModSeq="KGAGDGS[+80.0]DEEVDGKADGAEAKPAE";
			retentionTime=2607.104f;
			precursorCharge=4;
		} else if (false) {
			peptideModSeq="NGHDGDTHQEDDGEKS[+80.0]D";
			retentionTime=1495.3553f;
			precursorCharge=3;
		} else if (false) {
			peptideModSeq="NTPSQHSHSIQHS[+80.0]PER";
			retentionTime=1256.3296f;
			precursorCharge=4;
		} else if (true) {
			// IRS1
			peptideModSeq="KGS[+80.0]GDYMPMSPK";
			retentionTime=2949.1633f;
			precursorCharge=2;
		} else if (false) {
			// IRS1
			peptideModSeq="KGSGDYMPMS[+80.0]PK";
			retentionTime=2949.1633f;
			precursorCharge=2;
		} else if (true) {
			peptideModSeq="NTPSQHSHSIQHS[+80.0]PER";
			retentionTime=1256.3296f;
			precursorCharge=3;
		} else if (true) {
			peptideModSeq="NTPS[+80.0]QHSHSIQHSPER";
			retentionTime=1256.3296f;
			precursorCharge=3;
		} else if (false) {
			peptideModSeq="RPMEEDGEEKS[+80.0]PSK";
			retentionTime=1434.3873f;
			precursorCharge=3;
		} else if (false) {
			peptideModSeq="AVT[+80.0]PVPTKTEEVSNLK";
			retentionTime=3658.4482f;
			precursorCharge=3;
		} else if (false) {
			peptideModSeq="IDDRDS[+80.0]DEEGASDR";
			retentionTime=1606f;
			precursorCharge=2;
		} else if (true) {
			peptideModSeq="RAGDLLEDS[+80.0]PKRPK";
			retentionTime=36*60f;
			precursorCharge=3;
		}

		LibraryEntry libentry=library.getEntries(peptideModSeq, precursorCharge, false).get(0);
		
		libentry=libentry.updateRetentionTime(retentionTime);
		double precursorMz=parameters.getAAConstants().getChargedMass(peptideModSeq, precursorCharge);
		
		ArrayList<Stripe> stripes=stripefile.getStripes(precursorMz, 0, Float.MAX_VALUE, false);
		System.out.println(precursorMz+", "+stripes.size());
		ArrayList<String> permutations=PhosphoPermuter.getPermutations(peptideModSeq, PeptideModification.phosphorylation, parameters.getAAConstants());
		PhosphoLocalizationData actuallyPhosphoData=localizer.extractPhosphoFormsFromStripes(peptideModSeq, precursorMz, precursorCharge, permutations, retentionTime, stripes, true);

		System.out.println("Just off of localization ions");
		ArrayList<String> keys=new ArrayList<String>(actuallyPhosphoData.getPassingForms().keySet());
		for (String sequenceKey : keys) {

			XYPoint point=actuallyPhosphoData.getLocalizationScores().get(sequenceKey);
			float rt=(float)point.x;
			float localizationScore=(float)point.y;
			
			TransitionRefinementData data=actuallyPhosphoData.getPassingForms().get(sequenceKey);
			System.out.println(sequenceKey+"\t"+data.getApexRT()+"\t"+rt+"\t"+localizationScore);
			
			HashMap<String, HashMap<FragmentIon, XYTrace>> uniqueFragmentIons=actuallyPhosphoData.getUniqueFragmentIons();
			HashMap<String, HashMap<FragmentIon, XYTrace>> otherFragmentIons=actuallyPhosphoData.getOtherFragmentIons();
			HashMap<FragmentIon, XYTrace> uniqueFragments=uniqueFragmentIons.get(sequenceKey);
			HashMap<FragmentIon, XYTrace> otherFragments=new HashMap<FragmentIon, XYTrace>(otherFragmentIons.get(sequenceKey));

			HashMap<FragmentIon, XYTrace> allFragments=new HashMap<FragmentIon, XYTrace>();
			allFragments.putAll(uniqueFragments);
			allFragments.putAll(otherFragments);
			ArrayList<XYTrace> uniqueFragmentsList=new ArrayList<XYTrace>(allFragments.values());
			XYTraceInterface[] fragmentTraces=uniqueFragmentsList.toArray(new XYTrace[uniqueFragmentsList.size()]);

			Charter.launchChart(sequenceKey+" Retention Time (min)", "Intensity", false, new Dimension(1000, 400), fragmentTraces);
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
		
		HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> allVsUniqueList=actuallyPhosphoData.getScoreTraces();
		ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
		for (Entry<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> entry : allVsUniqueList.entrySet()) {
			String seq=entry.getKey();
			if (actuallyPhosphoData.getPassingForms().containsKey(seq)) {
				Pair<TFloatFloatHashMap, TFloatFloatHashMap> pair=entry.getValue();
				Color color=RandomGenerator.randomColor(seq.hashCode()*16807);
				//traces.add(new XYTrace(pair.x, GraphType.line, "ALL_"+seq, color, 5.0f));
				traces.add(new XYTrace(pair.y, GraphType.line, "UNI_"+seq, color, 3.0f));
			}
		}
		//traces.add(new XYTrace(primary, GraphType.boldline, "primary"));
		
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

		System.out.println("Based on all ions");
		ArrayList<LibraryEntry> entries=new ArrayList<>();
		entries.add(libentry);
		BlockingQueue<PeptideScoringResult> resultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
		BlockingQueue<ModificationLocalizationData> localizationQueue=new LinkedBlockingQueue<ModificationLocalizationData>();
		CASiLOneScoringTask task=new CASiLOneScoringTask(scorer, entries, stripes, dutyCycle, precursors, localizer, resultsQueue, localizationQueue, parameters);
		task.call();

		int index=0;
		while (!resultsQueue.isEmpty()) {
			if (!resultsQueue.isEmpty()) {
				PeptideScoringResult result=resultsQueue.take();
				ArrayList<Pair<ScoredObject<Stripe>, float[]>> data=result.getGoodStripes();
				index++;
				for (Pair<ScoredObject<Stripe>, float[]> pair : data) {
					System.out.println(index+") "+result.getEntry().getPeptideModSeq()+"\t"+pair.x.x+"\t("+((pair.x.y.getScanStartTime())/60f)+" minutes)");
				}
			} else {
				Thread.sleep(10);
			}
		}
		
		while(!localizationQueue.isEmpty()) {
			if (!localizationQueue.isEmpty()) {
				ModificationLocalizationData data=localizationQueue.take();
				System.out.println(data.getLocalizationPeptideModSeq().getPeptideAnnotation()+" ("+data.isSiteSpecific()+") --> "+data.getLocalizationScore()+"\t"+data.getLocalizingIntensity()+"\t"+data.getTotalIntensity());
			} else {
				Thread.sleep(10);
			}
		}
	}
}
