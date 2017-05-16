package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class PhosphoLocalizationScoringTask extends AbstractLibraryScoringTask {
	private final PhosphoLocalizer localizer;
	
	public PhosphoLocalizationScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, PrecursorScanMap precursors, PhosphoLocalizer localizer, BlockingQueue<PeptideScoringResult> resultsQueue,
			SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		this.localizer=localizer;
	}

	@Override
	protected Nothing process() {
		EncyclopediaScorer eScorer=(EncyclopediaScorer)scorer;
		for (LibraryEntry entry : super.entries) {
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq(), parameters.getAAConstants());
			AuxillaryPSMScorer auxScorer=eScorer.getAuxScorer().getEntryOptimizedScorer(entry);
			
			ArrayList<String> permutations=PhosphoPermuter.getPermutations(entry.getPeptideModSeq(), parameters.getAAConstants());
			PhosphoLocalizationData phosphoData=localizer.extractPhosphoFormsFromStripes(entry.getPeptideModSeq(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), permutations, entry.getRetentionTime(), super.stripes);
			
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
			
			for (String sequenceKey : keys) {
				String peptideModSeq=sequenceKey.replaceAll("\\(", "").replaceAll("\\)", "");
				
				XYPoint point=phosphoData.getLocalizationScores().get(sequenceKey);
				float rt=(float)point.x;
				float score=(float)point.y;
				
				FragmentationModel model=new FragmentationModel(peptideModSeq, parameters.getAAConstants());
				AnnotatedLibraryEntry unitEntry=model.getUnitSpectrum(entry.getSource(), entry.getAccessions(), entry.getPrecursorCharge(), rt, parameters, entry.isDecoy());

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
				
				float[] auxScoreArray=auxScorer.score(unitEntry, bestStripe, predictedIsotopeDistribution, precursors);
				float evalue=-1f;
				
				PeptideScoringResult result=new PeptideScoringResult(unitEntry);
				result.addStripe(score, General.concatenate(auxScoreArray, evalue), bestStripe);
				resultsQueue.add(result);
			}
		}
		return Nothing.NOTHING;
	}

	
}
