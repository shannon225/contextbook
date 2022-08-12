package edu.washington.gs.maccoss.encyclopedia.filereaders.spectrumprocessors;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map.Entry;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MSMSBlock;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SpectrumUtils;

public class WindowDownsampler implements Runnable {
	private final MassTolerance tolerance;
	private final ArrayList<Range> targetRanges;
	private final BlockingQueue<MSMSBlock> inputQueue;
	private final BlockingQueue<MSMSBlock> outputQueue;

	private Throwable error;

	public WindowDownsampler(ArrayList<Range> targetRanges, MassTolerance tolerance, BlockingQueue<MSMSBlock> inputQueue, BlockingQueue<MSMSBlock> outputQueue) {
		this.targetRanges=targetRanges;
		this.tolerance=tolerance;
		this.inputQueue=inputQueue;
		this.outputQueue=outputQueue;
	}

	public void run() {
		HashMap<Range, Range> targetRangeByActualRangeMap=new HashMap<>();
		HashSet<Range> allRanges=new HashSet<Range>();
		
		Range cycleStart=null;
		int cycleLength=0;

		try {
			LinkedList<FragmentScan> currentCycle=new LinkedList<FragmentScan>();
			while (true) {
				MSMSBlock block=inputQueue.take();
				if (MSMSBlock.POISON_BLOCK==block) {
					outputQueue.put(MSMSBlock.POISON_BLOCK);
					// needs to join up here if we start using multiple threads
					break;
				}
				
				// scan ahead to set up
				if (cycleStart==null) {
					STARTUP: for (FragmentScan stripe : block.getFragmentScans()) {
						Range thisRange=stripe.getRange();
						Range truncatedRange=new Range((int)thisRange.getStart(), (int)thisRange.getStop()); // to deal with rounding errors (works out to 1600 m/z)
						
						if (cycleStart==null) {
							cycleStart=truncatedRange;
						} else if (allRanges.contains(truncatedRange)) {
							break STARTUP;
						}

						cycleLength++;
						TARGETSCAN: for (Range target : targetRanges) {
							if (target.contains(truncatedRange)) {
								targetRangeByActualRangeMap.put(truncatedRange, target);
								break TARGETSCAN;
							}
						}
					}
					
				}
				
				ArrayList<FragmentScan> downsampledStripes=new ArrayList<FragmentScan>();
				BLOCK: for (FragmentScan stripe : block.getFragmentScans()) {
					currentCycle.add(stripe);
					if (currentCycle.size()>cycleLength) {
						currentCycle.removeFirst();
					} else if (currentCycle.size()<cycleLength) {
						continue BLOCK;
					}
					
					// process current cycle
					if (currentCycle.getFirst().getRange().equals(cycleStart)) {
						HashMap<Range, ArrayList<FragmentScan>> scansByTargetRange=new HashMap<>();
						for (FragmentScan thisStripe : currentCycle) {
							Range thisRange=thisStripe.getRange();
							Range truncatedRange=new Range((int)thisRange.getStart(), (int)thisRange.getStop()); // to deal with rounding errors (works out to 1600 m/z)
							
							Range targetRange=targetRangeByActualRangeMap.get(truncatedRange);
							
							ArrayList<FragmentScan> list=scansByTargetRange.get(targetRange);
							if (list==null) {
								list=new ArrayList<>();
								scansByTargetRange.put(targetRange, list);
							}
							list.add(thisStripe);
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
							downsampledStripes.add(newScan);
						}
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

	public boolean hadError() {
		return null != error;
	}

	public Throwable getError() {
		return error;
	}
}
