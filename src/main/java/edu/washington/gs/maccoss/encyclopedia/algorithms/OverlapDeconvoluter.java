package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlBlock;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import gnu.trove.list.array.TFloatArrayList;

public class OverlapDeconvoluter implements Runnable {
	private final MassTolerance tolerance;
	private final BlockingQueue<MzmlBlock> inputQueue;
	private final BlockingQueue<MzmlBlock> outputQueue;
	private final HashMap<Range, TFloatArrayList> retentionTimesByStripe=new HashMap<Range, TFloatArrayList>();

	public OverlapDeconvoluter(MassTolerance tolerance, BlockingQueue<MzmlBlock> inputQueue, BlockingQueue<MzmlBlock> outputQueue) {
		this.tolerance=tolerance;
		this.inputQueue=inputQueue;
		this.outputQueue=outputQueue;
	}
	
	public HashMap<Range, TFloatArrayList> getRetentionTimesByStripe() {
		return retentionTimesByStripe;
	}

	public void run() {
		Range cycleStart=null;
		int cycleLength=-1;
		int cycleLengthP2=-1;
		int cycleCenter=-1;
		
		try {
			LinkedList<Stripe> currentCycle=new LinkedList<Stripe>();
			while (true) {
				MzmlBlock block=inputQueue.take();
				if (MzmlBlock.POISON_BLOCK==block) {
					outputQueue.put(MzmlBlock.POISON_BLOCK);
					// needs to join up here if we start using multiple threads
					break;
				}
				if (cycleStart==null) {
					cycleStart=block.getStripes().get(0).getRange();
				}
				ArrayList<Stripe> deconvolutedStripes=new ArrayList<Stripe>();
				
				for (Stripe stripe : block.getStripes()) {
					if (cycleLength==-1) {
						if (currentCycle.size()>0&&cycleStart.equals(stripe.getRange())) {
							cycleLength=currentCycle.size();
							cycleLengthP2=cycleLength+2;
							cycleCenter=(int)Math.ceil(cycleLength/2.0f);
						}
						currentCycle.add(stripe);
						continue;
					}

					currentCycle.add(stripe);
					if (currentCycle.size()>cycleLengthP2) {
						currentCycle.removeFirst();
					}
					
					// START THREADABLE SECTION (should I need to)
					Stripe earlyLow=null;
					Stripe earlyHigh=null;
					Stripe center=null;
					Stripe lateLow=null;
					Stripe lateHigh=null;
					int count=0;
					for (Stripe s : currentCycle) {
						if (count==0) {
							earlyLow=s;
						} else if (count==1) {
							earlyHigh=s;
						} else if (count==cycleCenter) {
							center=s;
						} else if (count==cycleLengthP2-2) {
							lateLow=s;
						} else if (count==cycleLengthP2-1) {
							lateHigh=s;
						}
						count++;
					}
					
					if (!earlyLow.getRange().contains(center.getRange().getStart())) {
						earlyLow=null;
					}
					if (!lateLow.getRange().contains(center.getRange().getStart())) {
						lateLow=null;
					}
					
					if (!earlyHigh.getRange().contains(center.getRange().getStop())) {
						earlyHigh=null;
					}
					if (!lateHigh.getRange().contains(center.getRange().getStop())) {
						lateHigh=null;
					}
					
					Pair<Stripe, Stripe> deconvoluted=deconvolute(earlyLow, earlyHigh, center, lateLow, lateHigh, tolerance);
					// END THREADABLE SECTION
					
					deconvolutedStripes.add(deconvoluted.x);
					deconvolutedStripes.add(deconvoluted.y);
					
					addRetentionTime(deconvoluted.x);
					addRetentionTime(deconvoluted.y);
					
				}
				outputQueue.put(new MzmlBlock(block.getPrecursors(), deconvolutedStripes));
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("DIA writing interrupted!");
			Logger.errorException(ie);
		}
	}
	public void addRetentionTime(Stripe thisStripe) {
		Range range=thisStripe.getRange();
		TFloatArrayList stripeRTs=retentionTimesByStripe.get(range);
		if (stripeRTs==null) {
			stripeRTs=new TFloatArrayList();
			retentionTimesByStripe.put(range, stripeRTs);
		}
		stripeRTs.add(thisStripe.getScanStartTime());
	}
	
	public static Pair<Stripe, Stripe> deconvolute(Stripe earlyLow, Stripe earlyHigh, Stripe center, Stripe lateLow, Stripe lateHigh, MassTolerance tolerance) {
		float[] intensities=center.getIntensityArray();
		double[] masses=center.getMassArray();
		
		Range lowerRange=new Range(center.getRange().getStart(), center.getRange().getMiddle());
		Range upperRange=new Range(center.getRange().getMiddle(), center.getRange().getStop());
		
		ArrayList<Peak> lowerPeaks=new ArrayList<Peak>();
		ArrayList<Peak> upperPeaks=new ArrayList<Peak>();
		
		for (int i=0; i<masses.length; i++) {
			float earlyLowIntensity=getIntensity(tolerance, masses[i], earlyLow);
			float earlyHighIntensity=getIntensity(tolerance, masses[i], earlyHigh);
			float lateLowIntensity=getIntensity(tolerance, masses[i], lateLow);
			float lateHighIntensity=getIntensity(tolerance, masses[i], lateHigh);
			
			float totalLow=earlyLowIntensity+lateLowIntensity;
			float totalHigh=earlyHighIntensity+lateHighIntensity;
			float total=totalLow+totalHigh;
			
			if (total>0.0f) {
				float fractionLow=totalLow/total;
				if (fractionLow>0.0f) {
					float intensity=intensities[i]*fractionLow;
					lowerPeaks.add(new Peak(masses[i], intensity));
				}
				
				float fractionHigh=totalHigh/total;
				if (fractionHigh>0.0f) {
					float intensity=intensities[i]*fractionHigh;
					upperPeaks.add(new Peak(masses[i], intensity));
				}
			}
		}
		
		Stripe lowerStripe=getDeconvolutedStripe(center, lowerRange, lowerPeaks, false);
		Stripe upperStripe=getDeconvolutedStripe(center, upperRange, upperPeaks, true);
		
		Pair<Stripe, Stripe> deconvoluted=new Pair<Stripe, Stripe>(lowerStripe, upperStripe);
		return deconvoluted;
	}

	private static Stripe getDeconvolutedStripe(Stripe center, Range lowerRange, ArrayList<Peak> lowerPeaks, boolean useNegativeScanNumber) {
		Pair<double[], float[]> arrays=Peak.toArrays(lowerPeaks);
		int scanNumber=useNegativeScanNumber?(-center.getSpectrumIndex()):center.getSpectrumIndex();
		
		Stripe lowerStripe=new Stripe(center.getSpectrumName(), center.getPrecursorName(), scanNumber, center.getScanStartTime(), lowerRange.getStart(), lowerRange.getStop(), arrays.x, arrays.y);
		return lowerStripe;
	}

	private static float getIntensity(MassTolerance tolerance, double mass, Stripe stripe) {
		// if we're at a boundary, return minimum value. This means if the peak
		// is not in the other stripe, we get it. Otherwise, they essentially
		// get the intensity.
		if (stripe==null) return Float.MIN_VALUE;
		
		int[] indicies=tolerance.getIndicies(stripe.getMassArray(), mass);
		float intensity=0.0f;
		for (int j=0; j<indicies.length; j++) {
			intensity+=stripe.getIntensityArray()[indicies[j]];
		}
		return intensity;
	}
}
