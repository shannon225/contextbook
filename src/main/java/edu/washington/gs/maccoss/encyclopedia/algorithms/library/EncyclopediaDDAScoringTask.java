package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.Scribe;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class EncyclopediaDDAScoringTask extends AbstractLibraryScoringTask {
	// FIXME shouldn't be static (long term)
	private static final HashMap<String, float[]> isotopeDistributions=new HashMap<String, float[]>();
	EncyclopediaScorer scorerFunction;
	
	public EncyclopediaDDAScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, PrecursorScanMap precursors, BlockingQueue<AbstractScoringResult> resultsQueue, SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		scorerFunction=(EncyclopediaScorer)scorer;
	}

	@Override
	protected Nothing process() {
		for (FragmentScan msms : super.stripes) {
//			try {System.out.println("SLEEPING"); Thread.sleep(10000);}
//			catch (Exception e) {}
			ArrayList<ScoredIndex> goodHits=new ArrayList<ScoredIndex>();
			TFloatFloatHashMap map=new TFloatFloatHashMap();
			for (int i=0; i<super.entries.size(); i++) {
				LibraryEntry entry=super.entries.get(i);
				boolean match=parameters.getPrecursorTolerance().equals(entry.getPrecursorMZ(), msms.getPrecursorMZ());
				for (int j = 0; j < Scribe.NUMBER_OF_ISOTOPES_ABOVE_MONOISOTOPIC; j++) {
					double target = entry.getPrecursorMZ()+(j+1)*MassConstants.neutronMass/entry.getPrecursorCharge();
					match=match||parameters.getPrecursorTolerance().equals(target, msms.getPrecursorMZ());	
				}
				if (match) {
					float score=scorerFunction.score(entry, msms);
					goodHits.add(new ScoredIndex(score, i));
					map.put(i, score);
				}
			}
			if (map.size()==0) continue;
			
			EValueCalculator calculator=new EValueCalculator(map, 0f, 0.5f);
			int index=Math.round(calculator.getMaxRT());
			float score=calculator.getMaxRawScore();
			float evalue=calculator.getNegLnEValue();
			if (Float.isNaN(evalue)) {
				evalue=-1.0f;
			}
			
			LibraryEntry entry=super.entries.get(index);
			float[] predictedIsotopeDistribution=getIsotopeDistribution(entry);
			float[] auxScoreArray=scorerFunction.auxScore(entry, msms, predictedIsotopeDistribution, precursors);

			AbstractScoringResult result=new PeptideScoringResult(entry);
			result.addStripe(score, General.concatenate(auxScoreArray, evalue), msms);
			resultsQueue.add(result);
		}
		return Nothing.NOTHING;
	}

	public float[] getIsotopeDistribution(LibraryEntry entry) {
		float[] predictedIsotopeDistribution=isotopeDistributions.get(entry.getPeptideModSeq());
		if (predictedIsotopeDistribution==null) {
			predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq(), parameters.getAAConstants());
			isotopeDistributions.put(entry.getPeptideModSeq(), predictedIsotopeDistribution);
		}
		return predictedIsotopeDistribution;
	}

	
}
