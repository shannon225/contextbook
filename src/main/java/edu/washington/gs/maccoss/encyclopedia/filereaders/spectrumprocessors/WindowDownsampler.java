package edu.washington.gs.maccoss.encyclopedia.filereaders.spectrumprocessors;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MSMSBlock;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.WindowData;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SpectrumUtils;
import gnu.trove.list.array.TFloatArrayList;

public class WindowDownsampler implements SpectrumProcessor {
	private final MassTolerance tolerance;
	private final List<Range> downsampledRangeList;
	private HashMap<Range, Range> targetRangeByActualRangeMap;
	private BlockingQueue<MSMSBlock> inputQueue;
	private BlockingQueue<MSMSBlock> outputQueue;
	private final HashMap<Range, TFloatArrayList> retentionTimesByStripe=new HashMap<Range, TFloatArrayList>();
	private final HashMap<Range, TFloatArrayList> ionInjectionTimesByStripe=new HashMap<Range, TFloatArrayList>();
	private final HashMap<Range, TFloatArrayList> truncatedRetentionTimesByStripe=new HashMap<Range, TFloatArrayList>();
	private final HashMap<Range, TFloatArrayList> truncatedIonInjectionTimesByStripe=new HashMap<Range, TFloatArrayList>();

	private Throwable error;

	public WindowDownsampler(List<Range> downsampledRangeList, MassTolerance tolerance) {
		this.downsampledRangeList=downsampledRangeList;
		this.tolerance=tolerance;
	}
	
	@Override
	public void initialize(StripeFileInterface inputFile, SearchParameters params, BlockingQueue<MSMSBlock> inputQueue, BlockingQueue<MSMSBlock> outputQueue) {
		Map<Range, WindowData> thisRangeMap=inputFile.getRanges();
		targetRangeByActualRangeMap=mapRanges(thisRangeMap, downsampledRangeList);
		this.inputQueue=inputQueue;
		this.outputQueue=outputQueue;
	}
	
	@Override
	public HashMap<Range, TFloatArrayList> getIonInjectionTimesByStripe() {
		return ionInjectionTimesByStripe;
	}
	@Override
	public HashMap<Range, TFloatArrayList> getRetentionTimesByStripe() {
		return retentionTimesByStripe;
	}
	
	protected static HashMap<Range, Range> mapRanges(Map<Range, WindowData> thisRangeMap, List<Range> downsampledRangeList) {
		HashMap<Range, Range> targetRangeByActualRangeMap=new HashMap<Range, Range>();
		for (Range thisRange : thisRangeMap.keySet()) {			
			Range truncatedRange=new Range((int)thisRange.getStart(), (int)thisRange.getStop()); // to deal with rounding errors (works out to 1600 m/z)
			
			TARGETSCAN: for (Range target : downsampledRangeList) {
				if (target.contains(truncatedRange)) {
					targetRangeByActualRangeMap.put(truncatedRange, target);
					break TARGETSCAN;
				}
			}
		}
		return targetRangeByActualRangeMap;
	}

	@Override
	public void run() {
		try {
			HashSet<Range> currentWindowsSeen=new HashSet<Range>();
			LinkedList<FragmentScan> currentCycle=new LinkedList<FragmentScan>();
			while (true) {
				MSMSBlock block=inputQueue.take();
				if (MSMSBlock.POISON_BLOCK==block) {
					outputQueue.put(MSMSBlock.POISON_BLOCK);
					// needs to join up here if we start using multiple threads
					break;
				}
				
				ArrayList<FragmentScan> downsampledStripes=new ArrayList<FragmentScan>();
				BLOCK: for (FragmentScan stripe : block.getFragmentScans()) {
					// keep accumulating until we've seen a duplicate window
					if (!currentWindowsSeen.contains(stripe.getRange())) {
						currentWindowsSeen.add(stripe.getRange());
						currentCycle.add(stripe);
						continue BLOCK;
					}
					
					// if we see a duplicate window, then process the current one and clear for the next one
					LinkedList<FragmentScan> thisCycle=currentCycle;
					currentCycle=new LinkedList<>();
					currentWindowsSeen.clear();
					
					// process current cycle
					HashMap<Range, ArrayList<FragmentScan>> scansByTargetRange=new HashMap<>();
					for (FragmentScan thisStripe : thisCycle) {
						Range thisRange=thisStripe.getRange();
						Range truncatedRange=new Range((int)thisRange.getStart(), (int)thisRange.getStop()); // to deal with rounding errors (works out to 1600 m/z)
						
						Range targetRange=targetRangeByActualRangeMap.get(truncatedRange);
						
						if (targetRange!=null) {
							ArrayList<FragmentScan> list=scansByTargetRange.get(targetRange);
							if (list==null) {
								list=new ArrayList<>();
								scansByTargetRange.put(targetRange, list);
							}
							list.add(thisStripe);
						}
					}
					
					for (Entry<Range, ArrayList<FragmentScan>> entry : scansByTargetRange.entrySet()) {
						Range target=entry.getKey();
						ArrayList<FragmentScan> spectrumList = entry.getValue();
						Collections.sort(spectrumList);
						Spectrum downsampled=SpectrumUtils.accurateMergeSpectra(spectrumList, tolerance);
						
						// data taken from representative
						FragmentScan representative = spectrumList.get(0);
						String spectrumName="Merged_"+representative.getSpectrumName();
						String precursorName="Merged_"+representative.getPrecursorName();
						int spectrumIndex=representative.getSpectrumIndex();
						int fraction=representative.getFraction();
						byte charge=representative.getCharge();

						// data taken from overall range
						double isolationWindowLower=target.getStart();
						double isolationWindowUpper=target.getStop();
						double[] massArray=downsampled.getMassArray();
						float[] intensityArray=downsampled.getIntensityArray();
						
						// data aggregated-
						float scanStartTime=0.0f; // average
						float ionInjectionTime=0.0f; // sum
						for (FragmentScan scan : spectrumList) {
							scanStartTime+=scan.getScanStartTime();
							ionInjectionTime+=scan.getIonInjectionTime();
						}
						scanStartTime=scanStartTime/spectrumList.size();
						
						FragmentScan newScan=new FragmentScan(spectrumName, precursorName, spectrumIndex, scanStartTime, fraction, ionInjectionTime, isolationWindowLower, isolationWindowUpper, massArray, intensityArray, charge);
						addRetentionTime(newScan);
						downsampledStripes.add(newScan);
					}

				}
				outputQueue.put(new MSMSBlock(block.getPrecursorScans(), downsampledStripes));
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("DIA writing interrupted!");
			Logger.errorException(ie);
		} catch (Throwable t) {
			Logger.errorLine("Window downsampler failed!");
			Logger.errorException(t);

			this.error = t;
		}
	}

	@Override
	public boolean hadError() {
		return null != error;
	}

	@Override
	public Throwable getError() {
		return error;
	}
	
	private void addRetentionTime(FragmentScan thisStripe) {
		Range range=thisStripe.getRange();
		Range truncatedRange=new Range((int)range.getStart(), (int)range.getStop()); // to deal with rounding errors
		TFloatArrayList stripeRTs=truncatedRetentionTimesByStripe.get(truncatedRange);
		TFloatArrayList stripeIITs=truncatedIonInjectionTimesByStripe.get(truncatedRange);
		if (stripeRTs==null) {
			stripeRTs=new TFloatArrayList();
			truncatedRetentionTimesByStripe.put(truncatedRange, stripeRTs);
			retentionTimesByStripe.put(range, stripeRTs);
			stripeIITs=new TFloatArrayList();
			truncatedIonInjectionTimesByStripe.put(truncatedRange, stripeIITs);
			ionInjectionTimesByStripe.put(range, stripeIITs);
		}
		stripeRTs.add(thisStripe.getScanStartTime());
		stripeIITs.add(thisStripe.getIonInjectionTime());
	}
}
