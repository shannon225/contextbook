package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.awt.Color;
import java.awt.Dimension;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
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
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class PhosphoLocalizerExample {

	@SuppressWarnings("unused")
	public static void main(String[] args) throws Exception {
		//StripeFile.OPEN_IN_PLACE=true;
		LibraryFile.OPEN_IN_PLACE=true;
		
		//File libraryFile=new File("/Users/searleb/Documents/school/localization_manuscript/VillenJ_Exactive_HumanPhosphoproteome.dlib");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/22jun2016_mcf7_phospho_1a.dia");

		File libraryFile=new File("/Users/searleb/Documents/school/localization_manuscript/VillenJ_Exactive_HumanPhosphoproteome.dlib");
		//File libraryFile=new File("/Users/searleb/Documents/projects/phosphopedia/VillenJ_Exactive_HumanPhosphoproteome.dlib");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/22jun2016_mcf7_phospho_1a.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/22jun2016_mcf7_phospho_2a.dia");
		File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/22jun2016_mcf7_phospho_2b.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/22jun2016_mcf7_phospho_2c.dia");

		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/22jun2016_mcf7_phospho_1b.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/22jun2016_mcf7_phospho_1c.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/22jun2016_mcf7_phospho_1f.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/22jun2016_mcf7_phospho_2c.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/22jun2016_mcf7_phospho_5c.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/22jun2016_mcf7_phospho_5c.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/elibs/22jun2016_mcf7_phospho_3a_160627233451.mzML");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/prms/20160718_FU_bcs_4a_PRM.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/prms/20160718_FU_bcs_4b_PRM.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/prms/20160718_FU_bcs_4c_PRM.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/phospho_repeats/final_data/hela_repeats/20170430_HeLa_phosp_DIA_B_01_170506220515.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/phospho_repeats/final_data/hela_repeats/20170430_HeLa_phosp_DIA_B_02_170507024206.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/phospho_repeats/final_data/hela_repeats/20170430_HeLa_phosp_DIA_B_03_170507071858.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/phospho_repeats/final_data/hela_repeats/20170430_HeLa_phosp_DIA_B_04.dia");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/hela_phospho/110515_bcs_hela_phospho_starved_20mz_500_900.dia");

		//File libraryFile=new File("/Users/searleb/Documents/phospho_localization/data/VillenJ_Exactive_HumanPhosphoproteome.elib");
		//File diaFile=new File("/Volumes/BriansSSD/phospho/hela_repeats_prism/thesaurus_recalibrated_20p/20170430_HeLa_phosp_DIA_B_04.dia");
		
		//File libraryFile=new File("/Users/searleb/Documents/phospho_localization/data/VillenJ_Exactive_HumanPhosphoproteome.elib");
		//File diaFile=new File("/Users/searleb/Documents/phospho_localization/data/hela/110515_bcs_hela_phospho_starved_20mz_500_900.dia");
		
		//diaFile=new File("/Users/searleb/Documents/phospho_localization/hela_repeats/dia/20170430_HeLa_phosp_DIA_B_01_170506220515.dia");
		//diaFile=new File("/Users/searleb/Documents/phospho_localization/hela_repeats/dia/20170430_HeLa_phosp_DIA_B_02_170507024206.dia");
		//diaFile=new File("/Users/searleb/Documents/phospho_localization/hela_repeats/dia/20170430_HeLa_phosp_DIA_B_03_170507071858.dia");
		//diaFile=new File("/Users/searleb/Documents/phospho_localization/hela_repeats/dia/20170430_HeLa_phosp_DIA_B_04.dia");
		//libraryFile=new File("/Users/searleb/Documents/phospho_localization/hela_repeats/dia/20170430_HeLa_phosp_DIA_B_01_170506220515.dia.elib");
		//libraryFile=new File("/Users/searleb/Documents/phospho_localization/hela_repeats/dia/20170430_HeLa_phosp_DIA_B_02_170507024206.dia.elib");
		//libraryFile=new File("/Users/searleb/Documents/phospho_localization/hela_repeats/dia/20170430_HeLa_phosp_DIA_B_03_170507071858.dia.elib");
		//libraryFile=new File("/Users/searleb/Documents/phospho_localization/hela_repeats/dia/20170430_HeLa_phosp_DIA_B_04.dia.elib");

		diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/rosenberger/chludwig_K150309_013_SW_0-Pmix_George_dilution_0_in_human.dia");
		libraryFile=new File("/Users/searleb/Documents/school/localization_manuscript/rosenberger/psgs_standard_consensus_filtered.dlib");
		//diaFile=new File("/Users/searleb/Documents/phospho_localization/0.6.2/20170430_HeLa_phosp_DIA_B_01_170506220515.dia");
		//libraryFile=new File("/Users/searleb/Documents/phospho_localization/0.6.2/20170430_HeLa_phosp_DIA_B_01_170506220515.dia.elib");
		diaFile=new File("/Users/searleb/Documents/phospho_localization/rosenberger/chludwig_K150309_013_SW_0-Pmix_George_dilution_0_in_human.dia");
		libraryFile=new File("/Users/searleb/Documents/phospho_localization/rosenberger/psgs_standard_consensus_filtered.dlib");
		
		LibraryFile library=new LibraryFile();
		library.openFile(libraryFile);
		
		HashMap<String, String> defaults=SearchParameterParser.getDefaultParameters();
		defaults.put("-localizationModification", "Phosphorylation");
		defaults.put("-ptol", "16.67");
		defaults.put("-ftol", "16.67");
		defaults.put("-lftol", "16.67");
		//defaults.put("-frag", "yonly");
		//defaults.put("-scoringBreadthType", "uncal20");
		defaults.put("-scoringBreadthType", "recal");
		

		SearchParameters parameters=SearchParameterParser.parseParameters(defaults);
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters, true);
		
		UnitBackgroundFrequencyCalculator unitbackground=new UnitBackgroundFrequencyCalculator(0.01f);
		BackgroundFrequencyInterface background=unitbackground;
		//background=BackgroundFrequencyCalculator.generateBackground(stripefile);
		
		float duration=stripefile.getGradientLength()/20.0f;
		
		PhosphoLocalizer localizer=new PhosphoLocalizer(stripefile, PeptideModification.phosphorylation, background, parameters);
		
		String peptideModSeq;
		float retentionTime;
		byte precursorCharge;
		if (true) {
			peptideModSeq="RLSISS[+79.966331]LNPSNALSR[+10.008269]";
			retentionTime=3927.18115234375f;
			precursorCharge=3;
		} else if (true) {
			peptideModSeq="AASAPAKES[+80.0]PR";
			retentionTime=18.2f*60f;
			precursorCharge=2;
		} else if (false) {
			// rosenberger
			peptideModSeq="KSST[+79.966331]APDEISTSTK[+8.014199]";
			retentionTime=17.6f*60f;
			precursorCharge=3;
		} else if (false) {
			peptideModSeq="S[+80.0]FSKEVEER";
			retentionTime=2225.5823f;
			precursorCharge=2;
		} else if (false) {
			// targeted
			peptideModSeq="IHLGS[+80.0]SPK";
			retentionTime=34.7f*60f;
			precursorCharge=2;
			
		} else if (false) {
			// targeted
			peptideModSeq="MGSS[+80.0]PLEVPKPR";
			retentionTime=56.95f*60f;
			precursorCharge=2;
			
		} else if (false) {
			// targeted
			peptideModSeq="LGS[+80.0]PKPER";
			retentionTime=1315f;
			precursorCharge=2;
			
		} else if (false) {
			// targeted
			peptideModSeq="VDS[+80.0]PSHGLVTSSLC[+57.0]IPSPAR";
			retentionTime=5010f;
			precursorCharge=3;
			
		} else if (false) {
			peptideModSeq="GIAPAS[+80.0]PMLGNASNPNKADIPER";
			retentionTime=4189.1591796875f;
			precursorCharge=3;
		} else if (true) {
			// repeat 2 sp|P83731|RL24_HUMAN
			peptideModSeq="AITGAS[+80.0]LADIMAK";
			retentionTime=5680.037109375f;
			precursorCharge=2;
		} else if (false) {
			// repeat 1 (2 mins apart)
			peptideModSeq="LSGGLGAGS[+80.0]C[+57.0]R";
			retentionTime=2360.37475585937f;
			precursorCharge=2;
		} else if (false) {
			// hela starved (1 min apart)
			peptideModSeq="FGTFGGLGS[+80.0]K";
			retentionTime=88*60f;
			precursorCharge=2;
		} else if (false) {
			// hela starved (1 min apart)
			peptideModSeq="FGTFGGLGS[+80.0]K";
			retentionTime=88*60f;
			precursorCharge=2;
		} else if (false) {
			// repeat 4
			peptideModSeq="ATAPQTQHVS[+80.0]PMR";
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
		} else if (false) {
			// IRS1
			peptideModSeq="KGS[+80.0]GDYMPMS[+80.0]PK";
			retentionTime=2949.1633f;
			precursorCharge=2;
		} else if (true) {
			// IRS1
			peptideModSeq="KGS[+80.0]GDYMPMSPK";
			retentionTime=2949.1633f;
			precursorCharge=2;
		} else if (true) {
			// IRS1
			peptideModSeq="KGSGDYMPMS[+80.0]PK";
			retentionTime=2949.1633f;
			precursorCharge=2;
		} else if (false) {
			peptideModSeq="NTPSQHSHSIQHS[+80.0]PER";
			retentionTime=1256.3296f;
			precursorCharge=3;
		} else if (false) {
			peptideModSeq="NTPS[+80.0]QHSHSIQHSPER";
			retentionTime=1256.3296f;
			precursorCharge=3;
		} else if (false) {
			// from 110515_bcs_hela_phospho_starved_20mz_500_900
			peptideModSeq="RPMEEDGEEKS[+80.0]PSK";
			retentionTime=1434.3873f;
			precursorCharge=3;
		} else if (true) {
			peptideModSeq="TDGFAEAIHS[+80.0]PQVAGVPR";
			retentionTime=4580.765625f;
			precursorCharge=3;
		} else if (true) {
			peptideModSeq="MDS[+80.0]DEDEKEGEEEKVAK";
			retentionTime=1842.24536132813f;
			precursorCharge=3;
		} else if (true) {
			peptideModSeq="RAGDLLEDS[+80.0]PKRPK";
			retentionTime=36*60f;
			precursorCharge=3;
		}

		System.out.println("Looking up "+peptideModSeq+", +"+precursorCharge+" ("+(retentionTime/60f)+")");
		LibraryEntry libentry=library.getEntries(peptideModSeq, precursorCharge, false).get(0);
		
		libentry=libentry.updateRetentionTime(retentionTime);
		double precursorMz=parameters.getAAConstants().getChargedMass(peptideModSeq, precursorCharge);
		
		ArrayList<Stripe> stripes=stripefile.getStripes(precursorMz, 0, Float.MAX_VALUE, false);
		System.out.println(precursorMz+", "+stripes.size());
		ArrayList<String> permutations=PhosphoPermuter.getPermutations(peptideModSeq, PeptideModification.phosphorylation, parameters.getAAConstants());
		PhosphoLocalizationData actuallyPhosphoData=localizer.extractPhosphoFormsFromStripes(peptideModSeq, precursorMz, precursorCharge, permutations, retentionTime, stripes, true);

		System.out.println("Just off of localization ions");
		ArrayList<String> keys=new ArrayList<String>(actuallyPhosphoData.getPassingForms().keySet());
		HashMap<String, Pair<String, Float>> bestKeys=new HashMap<>();
		for (String sequenceKey : keys) {
			String pms=actuallyPhosphoData.getPassingForms().get(sequenceKey).getPeptideModSeq();
			XYPoint point=actuallyPhosphoData.getLocalizationScores().get(sequenceKey);
			float localizationScore=(float)point.y;
			if (bestKeys.containsKey(pms)) {
				Pair<String, Float> scored=bestKeys.get(pms);
				if (scored.y<localizationScore) {
					bestKeys.put(pms, new Pair<String, Float>(sequenceKey, localizationScore));
				}
			} else {
				bestKeys.put(pms, new Pair<String, Float>(sequenceKey, localizationScore));
			}
		}
		for (Pair<String, Float> pair : bestKeys.values()) {
			String sequenceKey=pair.x;
			XYPoint point=actuallyPhosphoData.getLocalizationScores().get(sequenceKey);
			float rt=(float)point.x;
			float localizationScore=(float)point.y;
			
			TransitionRefinementData data=actuallyPhosphoData.getPassingForms().get(sequenceKey);
			
			System.out.println(sequenceKey+"\t"+data.getApexRT()+"\t"+rt+"\t"+localizationScore);
			FragmentIon[] ions=data.getFragmentMassArray();
			float[] correlations=data.getCorrelationArray();
			
			
			HashMap<String, HashMap<FragmentIon, XYTrace>> uniqueFragmentIons=actuallyPhosphoData.getUniqueFragmentIons();
			HashMap<String, HashMap<FragmentIon, XYTrace>> otherFragmentIons=actuallyPhosphoData.getOtherFragmentIons();
			HashMap<FragmentIon, XYTrace> uniqueFragments=uniqueFragmentIons.get(sequenceKey);
			HashMap<FragmentIon, XYTrace> otherFragments=new HashMap<FragmentIon, XYTrace>(otherFragmentIons.get(sequenceKey));

			HashMap<FragmentIon, XYTrace> allFragments=new HashMap<FragmentIon, XYTrace>();
			for (int i=0; i<correlations.length; i++) {
				
				//if (correlations[i]>=TransitionRefiner.identificationCorrelationThreshold) {
					XYTrace unique=uniqueFragments.get(ions[i]);
					if (unique!=null) allFragments.put(ions[i], unique);
					XYTrace other=otherFragments.get(ions[i]);
					if (other!=null) allFragments.put(ions[i], other);
				//}
			}
			//allFragments.putAll(uniqueFragments);
			//allFragments.putAll(otherFragments);
			ArrayList<XYTrace> uniqueFragmentsList=new ArrayList<XYTrace>(allFragments.values());
			XYTraceInterface[] fragmentTraces=uniqueFragmentsList.toArray(new XYTrace[uniqueFragmentsList.size()]);

			Charter.launchChart(sequenceKey+" Retention Time (min)", "Intensity", false, new Dimension(500, 300), fragmentTraces);
		}

		EncyclopediaOneScorer scorer=new EncyclopediaOneScorer(parameters, unitbackground);
		FragmentationModel model=PeptideUtils.getPeptideModel(libentry.getPeptideModSeq(), parameters.getAAConstants());
		FragmentIon[] ions=model.getPrimaryIonObjects(parameters.getFragType(), libentry.getPrecursorCharge(), false, true);
		TFloatFloatHashMap primary=new TFloatFloatHashMap();
		for (int i=0; i<stripes.size(); i++) {
			Stripe stripe=stripes.get(i);
			if (stripe.getScanStartTime()>retentionTime-duration&&stripe.getScanStartTime()<retentionTime+duration) {
				primary.put(stripe.getScanStartTime()/60f, scorer.score(libentry, stripe, ions));
			}
		}

		ArrayList<Spectrum> limitedPrecursors=new ArrayList<>();
		for (PrecursorScan stripe : stripefile.getPrecursors(libentry.getRetentionTime()-duration, libentry.getRetentionTime()+duration)) {
			limitedPrecursors.add(stripe);
		}
		Charter.launchChart("Retention Time (min)", "Intensity", false, new Dimension(500, 300), ChromatogramExtractor.extractPrecursorChromatograms(parameters.getPrecursorTolerance(), libentry.getPrecursorMZ(), libentry.getPrecursorCharge(), limitedPrecursors));
		
		HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> allVsUniqueList=actuallyPhosphoData.getScoreTraces();
		ArrayList<XYTrace> traces=new ArrayList<XYTrace>();

		for (String seq : allVsUniqueList.keySet()) {
		//for (Pair<String, Float> pp : bestKeys.values()) {
		//	String seq=pp.x;
			System.out.println("HERE: "+seq);
			
			//if (actuallyPhosphoData.getPassingForms().containsKey(seq)) {
				Pair<TFloatFloatHashMap, TFloatFloatHashMap> pair=allVsUniqueList.get(seq);
				Color color=RandomGenerator.randomColor(seq.hashCode()*16807);
				//traces.add(new XYTrace(pair.x, GraphType.line, "ALL_"+seq, color, 5.0f));
				traces.add(new XYTrace(pair.y, GraphType.line, "UNI_"+seq, color, 3.0f));
			//}
		}
		//traces.add(new XYTrace(primary, GraphType.boldline, "primary"));
		
		Charter.launchChart("Retention Time (min)", "Localization Score", true, new Dimension(700, 300), traces.toArray(new XYTrace[traces.size()]));

		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));
		Range range=null;
		for (Range thisRange : stripefile.getRanges().keySet()) {
			if (thisRange.contains(precursorMz)) {
				range=thisRange;
				break;
			}
		}
		float dutyCycle=stripefile.getRanges().get(range);

		//TransitionRefiner.DISPLAY_PLOTS=true;
		System.out.println("Based on all ions");
		ArrayList<LibraryEntry> entries=new ArrayList<>();
		entries.add(libentry);
		BlockingQueue<PeptideScoringResult> resultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
		BlockingQueue<ModificationLocalizationData> localizationQueue=new LinkedBlockingQueue<ModificationLocalizationData>();
		ThesaurusOneScoringTask task=new ThesaurusOneScoringTask(scorer, entries, stripes, dutyCycle, precursors, localizer, resultsQueue, localizationQueue, parameters);
		//CASiLOneScoringTask task=new CASiLOneScoringTask(scorer, entries, stripes, dutyCycle, precursors, localizer, resultsQueue, localizationQueue, parameters);
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
				System.out.println(data.getLocalizationPeptideModSeq().getPeptideAnnotation()+" ("+data.isSiteSpecific()+","+data.isLocalized()+") --> "+data.getLocalizationScore()+"\t"+data.getLocalizingIntensity()+"\t"+data.getTotalIntensity()+"\t"+FragmentIon.toArchiveString(data.getLocalizingIons()));
				for (FragmentIon ion : data.getLocalizingIons()) {
					System.out.println("\t"+ion);
				}
			} else {
				Thread.sleep(10);
			}
		}
		
		/*ArrayList<XYTrace> primaryScoreTraces=new ArrayList<>();
		for (Entry<AmbiguousPeptideModSeq, TFloatFloatHashMap> entry : task.scoreByRTMapByPeptideAnnotation.entrySet()) {
			TFloatFloatHashMap map=entry.getValue();
			TFloatFloatHashMap sparse=new TFloatFloatHashMap();
			map.forEachEntry(new TFloatFloatProcedure() {
				@Override
				public boolean execute(float a, float b) {
					if (b>7) {
						sparse.put(a/60f, b);
					} else {
						sparse.put(a/60f, 7);
					}
					return true;
				}
			});
			Color color=RandomGenerator.randomColor(entry.getKey().getPeptideAnnotation().hashCode()*16807);
			XYTrace trace=new XYTrace(sparse, GraphType.line, entry.getKey().getPeptideAnnotation(), color, 3.0f);
			Charter.launchChart("Retention Time (min)", "Thesaurus Score", false, new Dimension(700, 200), trace);
			primaryScoreTraces.add(trace);
		}
		Charter.launchChart("Retention Time (min)", "Thesaurus Score", true, new Dimension(700, 300), primaryScoreTraces.toArray(new XYTrace[primaryScoreTraces.size()]));*/
	}
}
