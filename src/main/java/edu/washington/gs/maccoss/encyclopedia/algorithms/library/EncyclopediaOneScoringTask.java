package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.set.hash.TIntHashSet;

public class EncyclopediaOneScoringTask extends AbstractLibraryScoringTask {
	private final float dutyCycle;
	private final Range precursorIsolationRange;
	
	public EncyclopediaOneScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, Range precursorIsolationRange, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<AbstractScoringResult> resultsQueue,
			SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		this.dutyCycle=dutyCycle;
		this.precursorIsolationRange=new Range(precursorIsolationRange.getStart(), precursorIsolationRange.getStop());
	}
	
	private static final int peaksKept=5;

	@Override
	protected Nothing process() {
		EncyclopediaScorer eScorer=(EncyclopediaScorer)scorer;
		int movingAverageLength=Math.round(parameters.getExpectedPeakWidth()/dutyCycle);
		for (LibraryEntry entry : super.entries) {
			if (parameters.getTopNTargetsUsed()>0) {
				entry=entry.trimToNPeaks(parameters.getTopNTargetsUsed(), parameters.getAAConstants());
			}
			
			AuxillaryPSMScorer auxScorer=eScorer.getAuxScorer().getEntryOptimizedScorer(entry);
			FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
			FragmentIon[] ions=model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), true);
			Optional<FragmentIon[]> modificationSpecificIons;
			if (parameters.isVerifyModificationIons()) {
				modificationSpecificIons=model.getModificationSpecificIonObjects(precursorIsolationRange, parameters.getFragType(), entry.getPrecursorCharge(), true);
			} else {
				modificationSpecificIons=Optional.empty();
			}
			
			ions=FragmentIon.getUniqueFragments(ions, parameters.getFragmentTolerance()); // ensure that all ions are unique within tolerance
			
			AbstractScoringResult result=new PeptideScoringResult(entry);
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq(), parameters.getAAConstants());
			
			float[] primary=new float[super.stripes.size()];
			for (int i=0; i<super.stripes.size(); i++) {
				FragmentScan stripe=super.stripes.get(i);
				primary[i]=eScorer.score(entry, stripe, ions);
				
				if (modificationSpecificIons.isPresent()) {
					// if modified signal represents less than 25% of the score then don't trust it
					float scoreFromModIons=eScorer.score(entry, stripe, modificationSpecificIons.get());
					if (scoreFromModIons/primary[i]<0.25f) {
						primary[i]=0.0f;
					}
				}
			}
			
			//float[] averagePrimary=gaussianCenteredAverage(primary, movingAverageLength);

			TFloatFloatHashMap map=new TFloatFloatHashMap();
			ArrayList<ScoredIndex> goodStripes=new ArrayList<ScoredIndex>();
			for (int i=0; i<primary.length; i++) {
				goodStripes.add(new ScoredIndex(primary[i], i));
				map.put(i, primary[i]);
			}
			Collections.sort(goodStripes);

			EValueCalculator calculator=new EValueCalculator(map, 0f, 0.5f);

			TIntHashSet takenScans=new TIntHashSet();
			int identifiedPeaks=0;
			for (int i=goodStripes.size()-1; i>=0; i--) {
				float score=goodStripes.get(i).x;
				int index=goodStripes.get(i).y;
				if (takenScans.contains(index)) {
					continue;
					
				} else {	
					FragmentScan stripe=super.stripes.get(index);
					float[] auxScoreArray=auxScorer.score(entry, stripe, predictedIsotopeDistribution, precursors);
					float evalue=calculator.getNegLnEValue(score);
					if (Float.isNaN(evalue)) {
						evalue=-1.0f;
					}

					float deltaPrecursorMass=auxScorer.getParentDeltaMassIndex()>=0?auxScoreArray[auxScorer.getParentDeltaMassIndex()]:0.0f;
					float deltaFragmentMass=auxScorer.getFragmentDeltaMassIndex()>=0?auxScoreArray[auxScorer.getFragmentDeltaMassIndex()]:0.0f;
					
					result.addStripe(score, General.concatenate(auxScoreArray, evalue), deltaPrecursorMass, deltaFragmentMass, stripe);
					
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
			
			resultsQueue.add(result);
		}
		return Nothing.NOTHING;
	}

	
}
