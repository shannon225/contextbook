package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map.Entry;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.set.hash.TIntHashSet;

public class PhosphoEncyclopediaOneScoringTask extends AbstractLibraryScoringTask {
	private final PhosphoLocalizer localizer;
	private final float dutyCycle;
	
	public PhosphoEncyclopediaOneScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, float dutyCycle, PrecursorScanMap precursors, 
			PhosphoLocalizer localizer, BlockingQueue<PeptideScoringResult> resultsQueue, SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		this.dutyCycle=dutyCycle;
		this.localizer=localizer;
	}
	
	private static final int peaksKept=1;

	@Override
	protected Nothing process() {
		EncyclopediaScorer eScorer=(EncyclopediaScorer)scorer;
		int movingAverageLength=Math.round(parameters.getExpectedPeakWidth()/dutyCycle);
		for (LibraryEntry entry : super.entries) {
			AuxillaryPSMScorer auxScorer=eScorer.getAuxScorer().getEntryOptimizedScorer(entry);
			FragmentationModel model=new FragmentationModel(entry.getPeptideModSeq(), parameters.getAAConstants());
			FragmentIon[] ions=model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge());
			
			PeptideScoringResult result=new PeptideScoringResult(entry);
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq(), parameters.getAAConstants());
			
			float[] primary=new float[super.stripes.size()];
			for (int i=0; i<super.stripes.size(); i++) {
				Stripe stripe=super.stripes.get(i);
				primary[i]=eScorer.score(entry, stripe, ions);
			}
			
			float[] averagePrimary=gaussianCenteredAverage(primary, movingAverageLength);

			TFloatFloatHashMap scoreByRTMap=new TFloatFloatHashMap();
			ArrayList<ScoredIndex> goodStripes=new ArrayList<ScoredIndex>();
			for (int i=0; i<averagePrimary.length; i++) {
				goodStripes.add(new ScoredIndex(primary[i], i));
				scoreByRTMap.put(super.stripes.get(i).getScanStartTime(), primary[i]);
			}
			Collections.sort(goodStripes);

			EValueCalculator calculator=new EValueCalculator(scoreByRTMap);

			TIntHashSet takenScans=new TIntHashSet();
			int identifiedPeaks=0;
			for (int i=goodStripes.size()-1; i>=0; i--) {
				int index=goodStripes.get(i).y;
				if (takenScans.contains(index)) {
					continue;
					
				} else {
					Stripe stripe=super.stripes.get(index);
					
					considerLocalizations(entry, stripe.getScanStartTime(), eScorer, predictedIsotopeDistribution, auxScorer, calculator);
					
					// block out a 40 scan window
					int lowerWindow=index-2*movingAverageLength;
					int upperWindow=index+2*movingAverageLength;
					for (int j=lowerWindow; j<=upperWindow; j++) {
						takenScans.add(j);
					}
					
					if (identifiedPeaks>peaksKept) {
						// keep N+1 peaks
						break;
					}
					identifiedPeaks++;
				}
			}
		}
		return Nothing.NOTHING;
	}


	private void considerLocalizations(LibraryEntry entry, float retentionTime, EncyclopediaScorer eScorer, 
			float[] predictedIsotopeDistribution, AuxillaryPSMScorer auxScorer, EValueCalculator calculator) {
		ArrayList<String> permutations=PhosphoPermuter.getPermutations(entry.getPeptideModSeq(), parameters.getAAConstants());
		PhosphoLocalizationData phosphoData=localizer.extractPhosphoFormsFromStripes(entry.getPeptideModSeq(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), permutations, retentionTime, super.stripes, false);
		
		ArrayList<String> keys=new ArrayList<String>(phosphoData.getPassingForms().keySet());
		if (keys.size()==0) {
			String bestKey=null;
			float bestScore=-Float.MAX_VALUE;
			for (Entry<String, XYPoint> mapping : phosphoData.getLocalizationScores().entrySet()) {
				if (mapping.getValue().y>bestScore) {
					bestScore=(float)mapping.getValue().y;
					bestKey=mapping.getKey();
				}
			}
			if (bestKey!=null) {
				keys.add(bestKey);
			}
		}
		
		// FIXME THINK ABOUT HOW THIS MIGHT PUSH PEPTIDES OFF THEIR DETECTION AREA
		
		for (String sequenceKey : keys) {
			String peptideModSeq=sequenceKey.replaceAll("\\(", "").replaceAll("\\)", "").replaceAll("<", "").replaceAll(">", "");
			
			XYPoint point=phosphoData.getLocalizationScores().get(sequenceKey);
			float rt=(float)point.x;
			float localizationScore=(float)point.y;
			
			Pair<FragmentationModel, LibraryEntry> localizedForm=entry.getEntryFromNewSequence(peptideModSeq, entry.getAccessions(), entry.isDecoy(), parameters);
			FragmentationModel localizedModel=localizedForm.x;
			LibraryEntry localizedEntry=localizedForm.y;
			FragmentIon[] ions=localizedModel.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge());
			
			//FragmentationModel model=new FragmentationModel(peptideModSeq, parameters.getAAConstants());
			//AnnotatedLibraryEntry unitEntry=model.getUnitSpectrum(entry.getSource(), entry.getAccessions(), entry.getPrecursorCharge(), rt, parameters, entry.isDecoy());

			Stripe bestStripe=null;
			float bestDeltaRT=Float.MAX_VALUE;
			for (int i=0; i<super.stripes.size(); i++) {
				Stripe stripe=super.stripes.get(i);
				float delta=Math.abs(stripe.getScanStartTime()-rt);
				if (delta<bestDeltaRT) {
					bestDeltaRT=delta;
					bestStripe=stripe;
				}
			}

			float[] auxScoreArray=auxScorer.score(localizedEntry, bestStripe, predictedIsotopeDistribution, precursors);

			float score=eScorer.score(localizedEntry, bestStripe, ions);
			float evalue=calculator.getNegLog10EValue(score);
			if (Float.isNaN(evalue)) {
				evalue=-1.0f;
			}
			
			PeptideScoringResult result=new PeptideScoringResult(localizedEntry);
			result.addStripe(score, General.concatenate(auxScoreArray, evalue, localizationScore), bestStripe);
			resultsQueue.add(result);
		}
	}
	
}
