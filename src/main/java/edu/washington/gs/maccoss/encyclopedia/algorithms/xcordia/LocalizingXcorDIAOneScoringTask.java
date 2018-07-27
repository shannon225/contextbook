package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.AmbiguousPeptideModSeq;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.CASiLOneScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.FragmentIonBlacklist;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoPermuter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusOneScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusOneScoringTask.LocalizableForm;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.VariantFastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntRange;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.ComparablePair;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.IndexedObject;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class LocalizingXcorDIAOneScoringTask extends AbstractLibraryScoringTask {
	private final float dutyCycle;
	private final BackgroundFrequencyInterface background;
	private final BlockingQueue<ModificationLocalizationData> localizationQueue;
	private final float minimumScore;
	private final int movingAverageLength;
	private final PeptideModification localizingModification;
	
	public LocalizingXcorDIAOneScoringTask(PSMScorer scorer, BackgroundFrequencyInterface background, ArrayList<LibraryEntry> entries, 
			ArrayList<Stripe> stripes, Range precursorIsolationRange, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue,
			BlockingQueue<ModificationLocalizationData> localizationQueue, SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		this.background=background;
		this.dutyCycle=dutyCycle;
		this.localizationQueue=localizationQueue;
		this.minimumScore=-Log.log10(parameters.getPercolatorThreshold());
		this.movingAverageLength=Math.round(parameters.getExpectedPeakWidth()/dutyCycle);
		if (parameters.getLocalizingModification().isPresent()) {
			this.localizingModification=parameters.getLocalizingModification().get();
		} else {
			this.localizingModification=PeptideModification.polymorphism; 
		}
		
		assert(scorer instanceof XCorDIAOneScorer);
	}

	@Override
	protected Nothing process() {
		// separate targets from decoys and process in batches
		ArrayList<XCorrLibraryEntry> targetBatch=new ArrayList<>();
		ArrayList<XCorrLibraryEntry> decoyBatch=new ArrayList<>();
		for (LibraryEntry entry : super.entries) {
			XCorrLibraryEntry xcordiaEntry = getXCorrEntry(entry, parameters);
			if (entry.isDecoy()) {
				decoyBatch.add(xcordiaEntry);
			} else {
				targetBatch.add(xcordiaEntry);
			}
		}
		processPeptide(targetBatch);
		processPeptide(decoyBatch);

		return Nothing.NOTHING;
	}
	
	private void processPeptide(ArrayList<XCorrLibraryEntry> seedEntries) {
		// score N peptides individually
		ArrayList<ScoredObject<IndexedObject<PeptidePrecursor>>> bestScoresByEntry=new ArrayList<>();
		HashMap<PeptidePrecursor, float[]> scoresByEntry=new HashMap<>();
		HashMap<PeptidePrecursor, float[]> isotopesByEntry=new HashMap<>();
		for (XCorrLibraryEntry xcordiaEntry : seedEntries) {
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(xcordiaEntry.getPeptideModSeq(), parameters.getAAConstants());
			isotopesByEntry.put(xcordiaEntry, predictedIsotopeDistribution);
			
			// perform initial scoring
			float[] primary=new float[super.stripes.size()];
			for (int i=0; i<super.stripes.size(); i++) {
				Stripe stripe=super.stripes.get(i);
				XCorrStripe xcordiaStripe;
				if (stripe instanceof XCorrStripe) {
					xcordiaStripe=(XCorrStripe)stripe;
				} else {
					xcordiaStripe=new XCorrStripe(stripe, parameters);
				}
				primary[i]=scorer.score(xcordiaEntry, xcordiaStripe, predictedIsotopeDistribution, precursors);
			}
			
			float[] averagePrimary=gaussianCenteredAverage(primary, movingAverageLength);
			
			// determine the N best peaks for this peptide
			LinkedList<ScoredIndex> keptIndicies=new LinkedList<>();
			for (int i=0; i<averagePrimary.length; i++) {
				if (keptIndicies.size()==0) {
					keptIndicies.add(new ScoredIndex(averagePrimary[i], i));
				} else if (averagePrimary[i]>keptIndicies.get(keptIndicies.size()-1).x) {
					int upperIndexRange=i+movingAverageLength;
					int lowerIndexRange=i-movingAverageLength;
					int count=0;
					for (ScoredIndex scoredIndex : keptIndicies) {
						if (averagePrimary[i]>scoredIndex.x&&(scoredIndex.y<lowerIndexRange||scoredIndex.y>upperIndexRange)) {
							break;
						}
						count++;
					}
					keptIndicies.add(count, new ScoredIndex(averagePrimary[i], i));
					if (keptIndicies.size()>seedEntries.size()) {
						keptIndicies.removeLast();
					}
				}
			}
			
			for (ScoredIndex scoredIndex : keptIndicies) {
				bestScoresByEntry.add(new ScoredObject<IndexedObject<PeptidePrecursor>>(scoredIndex.x, new IndexedObject<PeptidePrecursor>(scoredIndex.y, xcordiaEntry)));
			}
			scoresByEntry.put(xcordiaEntry, averagePrimary);
		}
		Collections.sort(bestScoresByEntry);
		
		// determine N total peaks worth scoring (at least one for each peptide)
		TIntArrayList pickedPeakIndicies=new TIntArrayList();
		HashSet<PeptidePrecursor> previouslyPicked=new HashSet<>();
		for (int i=bestScoresByEntry.size()-1; i>=0; i--) {
			ScoredObject<IndexedObject<PeptidePrecursor>> scoredPeptide=bestScoresByEntry.get(i);
			
			if (!previouslyPicked.contains(scoredPeptide.y.y)) {
				previouslyPicked.add(scoredPeptide.y.y);
				pickedPeakIndicies.add(scoredPeptide.y.x);
			}
			
			if (previouslyPicked.size()==seedEntries.size()) break;
		}
		
		// localize each of N peaks
		for (int index : pickedPeakIndicies.toArray()) {
			
			// check quant
		}
	}

	private static XCorrLibraryEntry getXCorrEntry(LibraryEntry entry, SearchParameters parameters) {
		XCorrLibraryEntry xcordiaEntry;
		if (entry instanceof XCorrLibraryEntry) {
			xcordiaEntry=(XCorrLibraryEntry)entry;
		} else {
			FastaPeptideEntry peptide=new FastaPeptideEntry(entry.getSource(), entry.getAccessions(), entry.getPeptideModSeq());
			xcordiaEntry=XCorrLibraryEntry.generateEntry(false, peptide, entry.getPrecursorCharge(), parameters);
		}
		xcordiaEntry.init();
		return xcordiaEntry;
	}
}
