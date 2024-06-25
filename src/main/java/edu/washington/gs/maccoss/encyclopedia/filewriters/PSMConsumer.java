package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ScoredPSM;
import edu.washington.gs.maccoss.encyclopedia.algorithms.SimplePSM;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SimplePeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class PSMConsumer implements PeptideScoringResultsConsumer {
	// only written to by one thread
	private final ArrayList<PSMInterface> savedResults=new ArrayList<PSMInterface>();
	private final BlockingQueue<AbstractScoringResult> resultsQueue;
	private final AminoAcidConstants aaConstants;

	public PSMConsumer(BlockingQueue<AbstractScoringResult> resultsQueue, AminoAcidConstants aaConstants) {
		this.resultsQueue=resultsQueue;
		this.aaConstants=aaConstants;
	}
	
	public ArrayList<PSMInterface> getSavedResults() {
		return new ArrayList<PSMInterface>(savedResults);
	}

	@Override
	public void close() {
	}

	@Override
	public int getNumberProcessed() {
		// may be out of date, but at least won't throw concurrency errors
		return savedResults.size();
	}

	@Override
	public BlockingQueue<AbstractScoringResult> getResultsQueue() {
		return resultsQueue;
	}

	@Override
	public void run() {
		try {
			while (true) {
				AbstractScoringResult result=resultsQueue.take();
				if (AbstractScoringResult.POISON_RESULT==result) break;
				
				if (!result.hasScoredResults()) continue;
				ScoredPSM scoredMSMS = result.getScoredMSMS();
				
				SimplePeptidePrecursor precursor=new SimplePeptidePrecursor(scoredMSMS.getLibraryEntry().getPeptideModSeq(), scoredMSMS.getLibraryEntry().getPrecursorCharge(), aaConstants);
				SimplePSM psm=new SimplePSM(precursor, scoredMSMS.getDeltaPrecursorMass(), scoredMSMS.getDeltaFragmentMass(), scoredMSMS.getRTData());
				savedResults.add(psm);
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("DIA processing interrupted!");
			Logger.errorException(ie);
		}

	}
}
