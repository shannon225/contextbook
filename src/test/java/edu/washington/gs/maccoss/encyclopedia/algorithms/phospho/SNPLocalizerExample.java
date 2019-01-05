package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.awt.Color;
import java.awt.Dimension;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
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
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class SNPLocalizerExample {
	
	public static void main2(String[] args) throws Exception {
		File diaFile=new File("/Volumes/DataBackup/23aug2017_hela_serum_timecourse/23aug2017_fingerprint_s_dia_002.mzML");

		HashMap<String, String> defaults=SearchParameterParser.getDefaultParameters();
		defaults.put("-localizationModification", "Phosphorylation");
		defaults.put("-ptol", "10");
		defaults.put("-ftol", "10");
		defaults.put("-lftol", "10");
		defaults.put("-scoringBreadthType", "uncal20");
		
		SearchParameters parameters=SearchParameterParser.parseParameters(defaults);
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		UnitBackgroundFrequencyCalculator unitBackgroundFrequencyCalculator=new UnitBackgroundFrequencyCalculator(0.01f);
		BackgroundFrequencyInterface background=unitBackgroundFrequencyCalculator;
		background=BackgroundFrequencyCalculator.generateBackground(stripefile);
		
		PhosphoLocalizer localizer=new PhosphoLocalizer(stripefile, PeptideModification.polymorphism, background, parameters);
		
		String peptideModSeq="SLDLDSIINAVR";
		float retentionTime=86.5f*60;
		byte precursorCharge=2;

		FragmentationModel model=PeptideUtils.getPeptideModel(peptideModSeq, parameters.getAAConstants());
		HashSet<String> accessions=new HashSet<>();
		accessions.add("KRT13");
		LibraryEntry libentry=model.getUnitSpectrum("23aug2017_fingerprint_p_dia_002", accessions, precursorCharge, retentionTime, parameters);
		
		libentry=libentry.updateRetentionTime(retentionTime);
		double precursorMz=parameters.getAAConstants().getChargedMass(peptideModSeq, precursorCharge);
		
		ArrayList<FragmentScan> stripes=stripefile.getStripes(precursorMz, 0, Float.MAX_VALUE, false);
		System.out.println(precursorMz+", "+stripes.size());
		ArrayList<String> permutations=PhosphoPermuter.getPermutations(peptideModSeq, PeptideModification.polymorphism, parameters.getAAConstants());
		permutations.add("SLDLDSIINAVR");
		permutations.add("SLDLDSIIDAVR"); 
		PhosphoLocalizationData actuallyPhosphoData=localizer.extractPhosphoFormsFromStripes(peptideModSeq, precursorMz, precursorCharge, permutations, retentionTime, stripes, true);

		System.out.println("Just off of localization ions");
		ArrayList<String> keys=new ArrayList<String>(actuallyPhosphoData.getLocalizationScores().keySet());
		for (String sequenceKey : keys) {

			XYPoint point=actuallyPhosphoData.getLocalizationScores().get(sequenceKey);
			float rt=(float)point.x;
			float localizationScore=(float)point.y;
			
			//TransitionRefinementData data=actuallyPhosphoData.getPassingForms().get(sequenceKey);
			//System.out.println(sequenceKey+"\t"+data.getApexRT()+"\t"+rt+"\t"+localizationScore);
			
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

		float duration=stripefile.getGradientLength()/20.0f;
		EncyclopediaOneScorer scorer=new EncyclopediaOneScorer(parameters, (UnitBackgroundFrequencyCalculator)unitBackgroundFrequencyCalculator);
		model=PeptideUtils.getPeptideModel(libentry.getPeptideModSeq(), parameters.getAAConstants());
		FragmentIon[] ions=model.getPrimaryIonObjects(parameters.getFragType(), libentry.getPrecursorCharge(), false, true);
		TFloatFloatHashMap primary=new TFloatFloatHashMap();
		for (int i=0; i<stripes.size(); i++) {
			FragmentScan stripe=stripes.get(i);
			if (stripe.getScanStartTime()>retentionTime-duration&&stripe.getScanStartTime()<retentionTime+duration) {
				primary.put(stripe.getScanStartTime()/60f, scorer.score(libentry, stripe, ions));
			}
		}
		
		HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> allVsUniqueList=actuallyPhosphoData.getScoreTraces();
		ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
		for (Entry<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> entry : allVsUniqueList.entrySet()) {
			String seq=entry.getKey();
			Pair<TFloatFloatHashMap, TFloatFloatHashMap> pair=entry.getValue();
			Color color=RandomGenerator.randomColor(seq.hashCode());
			//traces.add(new XYTrace(pair.x, GraphType.line, "ALL_"+seq, color, 5.0f));
			traces.add(new XYTrace(pair.y, GraphType.line, "UNI_"+seq, color, 3.0f));
		}
		//traces.add(new XYTrace(primary, GraphType.boldline, "primary"));
		
		Charter.launchChart("Retention Time (All Ions)", "Localization Score", true, new Dimension(1000, 400), traces.toArray(new XYTrace[traces.size()]));
	}

	public static void main(String[] args) throws Exception {
		//File diaFile=new File("/Volumes/BriansSSD/20150908_6BB2_DIA_01.dia ");
		//File diaFile=new File("/Users/searleb/Documents/school/xcordia_manuscript/demux/20141121_3_4_DIA_1.dia");
		File diaFile=new File("/Users/searleb/Documents/backup/xcordia_manuscript/xcordia_5p/20141121_3_4_DIA_1.dia");

		HashMap<String, String> defaults=SearchParameterParser.getDefaultParameters();
		defaults.put("-localizationModification", "Phosphorylation");
		defaults.put("-ptol", "16.67");
		defaults.put("-ftol", "16.67");
		defaults.put("-lftol", "16.67");
		defaults.put("-scoringBreadthType", "uncal20");
		
		SearchParameters parameters=SearchParameterParser.parseParameters(defaults);
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		UnitBackgroundFrequencyCalculator unitBackgroundFrequencyCalculator=new UnitBackgroundFrequencyCalculator(0.01f);
		BackgroundFrequencyInterface background=unitBackgroundFrequencyCalculator;
		//background=BackgroundFrequencyCalculator.generateBackground(stripefile);
		
		PhosphoLocalizer localizer=new PhosphoLocalizer(stripefile, PeptideModification.polymorphism, background, parameters);
		
		String peptideModSeq="GPGGAWAAEVISDAR";
		float retentionTime=43*60;
		byte precursorCharge=2;

		FragmentationModel model=PeptideUtils.getPeptideModel(peptideModSeq, parameters.getAAConstants());
		HashSet<String> accessions=new HashSet<>();
		accessions.add("SAA1_HUMAN");
		accessions.add("SAA2_HUMAN");
		LibraryEntry libentry=model.getUnitSpectrum("20150908_6BB2_DIA_01", accessions, precursorCharge, retentionTime, parameters);
		
		libentry=libentry.updateRetentionTime(retentionTime);
		double precursorMz=parameters.getAAConstants().getChargedMass(peptideModSeq, precursorCharge);
		
		ArrayList<FragmentScan> stripes=stripefile.getStripes(precursorMz, 0, Float.MAX_VALUE, false);
		System.out.println(precursorMz+", "+stripes.size());
		ArrayList<String> permutations=PhosphoPermuter.getPermutations(peptideModSeq, PeptideModification.polymorphism, parameters.getAAConstants());
		//permutations.add("GPGGVWAAEAISDAR"); // 728.8626 M/Z SAA1
		//permutations.add("GPGGAWAAEVISNAR"); // 728.3706 M/Z SAA2
		//permutations.add("GPGGAWAAEVISDAR"); // 728.8626 M/Z SAA1.5
		permutations.add("GPGGAWAAEVISDAR");
		permutations.add("GPGGAWAAEVISTAR");
		permutations.add("GPGGVWAAEAISDAR");
		permutations.add("GPGGVWAAEAISNAR");
		permutations.add("GPGGAWAAEVISNAR");
		permutations.add("GPGGAWAADVISNAR");

		PhosphoLocalizationData actuallyPhosphoData=localizer.extractPhosphoFormsFromStripes(peptideModSeq, precursorMz, precursorCharge, permutations, retentionTime, stripes, true);

		System.out.println("Just off of localization ions");
		ArrayList<String> keys=new ArrayList<String>(actuallyPhosphoData.getLocalizationScores().keySet());
		for (String sequenceKey : keys) {

			XYPoint point=actuallyPhosphoData.getLocalizationScores().get(sequenceKey);
			float rt=(float)point.x;
			float localizationScore=(float)point.y;
			
			TransitionRefinementData data=actuallyPhosphoData.getPassingForms().get(sequenceKey);
			if (data==null) {
				continue;
			}
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

			if ("GPGGVWAAEAISDAR".equals(data.getPeptideModSeq())||"GPGGAWAAEVISDAR".equals(data.getPeptideModSeq())) {
				Charter.launchChart(sequenceKey+" Retention Time (min)", "Intensity", false, new Dimension(450, 300), fragmentTraces);
			}
		}

//		float duration=stripefile.getGradientLength()/20.0f;
//		EncyclopediaOneScorer scorer=new EncyclopediaOneScorer(parameters, (UnitBackgroundFrequencyCalculator)unitBackgroundFrequencyCalculator);
//		model=PeptideUtils.getPeptideModel(libentry.getPeptideModSeq(), parameters.getAAConstants());
//		FragmentIon[] ions=model.getPrimaryIonObjects(parameters.getFragType(), libentry.getPrecursorCharge());
//		TFloatFloatHashMap primary=new TFloatFloatHashMap();
//		for (int i=0; i<stripes.size(); i++) {
//			Stripe stripe=stripes.get(i);
//			if (stripe.getScanStartTime()>retentionTime-duration&&stripe.getScanStartTime()<retentionTime+duration) {
//				primary.put(stripe.getScanStartTime()/60f, scorer.score(libentry, stripe, ions));
//			}
//		}
		
		HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> allVsUniqueList=actuallyPhosphoData.getScoreTraces();
		ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
		for (Entry<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> entry : allVsUniqueList.entrySet()) {
			String seq=entry.getKey();
			Pair<TFloatFloatHashMap, TFloatFloatHashMap> pair=entry.getValue();
			Color color=RandomGenerator.randomColor(seq.hashCode());
			//traces.add(new XYTrace(pair.x, GraphType.line, "ALL_"+seq, color, 5.0f));
			traces.add(new XYTrace(pair.y, GraphType.line, "UNI_"+seq, color, 3.0f));
		}
		//traces.add(new XYTrace(primary, GraphType.boldline, "primary"));
		
		Charter.launchChart("Retention Time (All Ions)", "Localization Score", true, new Dimension(1000, 400), traces.toArray(new XYTrace[traces.size()]));

	}
}