package edu.washington.gs.maccoss.encyclopedia.filereaders.spectrumprocessors;

import java.util.HashMap;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MSMSBlock;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import gnu.trove.list.array.TFloatArrayList;

public interface SpectrumProcessor extends Runnable {
	public void initialize(StripeFileInterface inputFile, SearchParameters params, BlockingQueue<MSMSBlock> inputQueue, BlockingQueue<MSMSBlock> outputQueue);

	public HashMap<Range, TFloatArrayList> getRetentionTimesByStripe();
	public HashMap<Range, TFloatArrayList> getIonInjectionTimesByStripe();
	
	public boolean hadError();
	public Throwable getError();
}
