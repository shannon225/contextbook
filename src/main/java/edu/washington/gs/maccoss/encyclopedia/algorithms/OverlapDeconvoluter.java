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
		boolean isEven=true;

		// to keep track of whether we're on the early block or the late block
		Range minimumRange=null;
		Range maximumRange=null;
		Range previousRange=new Range(Float.MAX_VALUE, Float.MAX_VALUE);
		boolean isOnEarlyBlock=false;

		try {
			LinkedList<Pair<Stripe, Boolean>> currentCycle=new LinkedList<Pair<Stripe, Boolean>>();
			while (true) {
				MzmlBlock block=inputQueue.take();
				if (MzmlBlock.POISON_BLOCK==block) {
					outputQueue.put(MzmlBlock.POISON_BLOCK);
					// needs to join up here if we start using multiple threads
					break;
				}
				
				// scan ahead to set up
				if (cycleStart==null) {
					STARTUP: for (Stripe stripe : block.getStripes()) {
						if (cycleStart==null) {
							cycleStart=block.getStripes().get(0).getRange();
							minimumRange=cycleStart;
							maximumRange=cycleStart;
							cycleLength++;
							continue STARTUP;
						}
						
						// FYI: this approach burns half of the first duty cycle
						// and half of the last one. We could try to recover
						// that data but it's probably not worth anything
						if (stripe.getRange().compareTo(minimumRange)<0) {
							minimumRange=stripe.getRange();
						} else if (stripe.getRange().compareTo(maximumRange)>0) {
							maximumRange=stripe.getRange();
						}
						
						cycleLength++;
						if (cycleStart.equals(stripe.getRange())) {
							cycleLengthP2=cycleLength+2;
							cycleCenter=(int)Math.ceil(cycleLength/2.0f);
							isEven=cycleLength%2==0;
							break STARTUP;
						}
						
					}
				}
				
				ArrayList<Stripe> deconvolutedStripes=new ArrayList<Stripe>();
				BLOCK: for (Stripe stripe : block.getStripes()) {
					if (stripe.getRange().compareTo(previousRange)<0) {
						// cycled back
						if (stripe.getRange().compareTo(minimumRange)==0) {
							// early block
							isOnEarlyBlock=true;
						} else {
							isOnEarlyBlock=false;
						}
					}

					currentCycle.add(new Pair<Stripe, Boolean>(stripe, isOnEarlyBlock));
					previousRange=stripe.getRange();
					if (currentCycle.size()>cycleLengthP2) {
						currentCycle.removeFirst();
					} else if (currentCycle.size()<cycleLengthP2) {
						continue BLOCK;
					}
					
					boolean firstIsOnEarlyBlock=true;
					boolean secondIsOnEarlyBlock=true;
					Stripe earlyLow=null;
					Stripe earlyHigh=null;
					Stripe center=null;
					Stripe centerPlusOne=null;
					Stripe lateLow=null;
					Stripe lateHigh=null;
					int count=0;
					for (Pair<Stripe, Boolean> p : currentCycle) {
						Stripe s=p.x;
						if (count==0) {
							earlyLow=s;
							firstIsOnEarlyBlock=p.y;
						} else if (count==1) {
							earlyHigh=s;
							secondIsOnEarlyBlock=p.y;
						} else if (count==cycleCenter) {
							center=s;
						} else if (count==cycleCenter+1) {
							centerPlusOne=s;
						} else if (count==cycleLengthP2-2) {
							lateLow=s;
						} else if (count==cycleLengthP2-1) {
							lateHigh=s;
						}
						count++;
					}

					// START THREADABLE SECTION (should I need to)
					ArrayList<Pair<Stripe, Stripe>> deconvoluted=new ArrayList<Pair<Stripe,Stripe>>();
					if (!isEven) {
						// odd logic for choosing the center is easy...
						deconvoluted.add(deconvolute(earlyLow, earlyHigh, center, lateLow, lateHigh, tolerance));
					} else {
						// but bin splitting logic for even bin numbers is a little crazy!
						if (earlyLow.getRange().compareTo(maximumRange)==0&&earlyHigh.getRange().compareTo(minimumRange)==0) {
							// spanning outside, so ignore
						} else if (earlyLow.getRange().compareTo(minimumRange)==0) {
							// center is left side of right block so C
							deconvoluted.add(deconvolute(earlyLow, earlyHigh, center, lateLow, lateHigh, tolerance));
						} else if (earlyHigh.getRange().compareTo(maximumRange)==0) {
							// center is right side of left block so C+1
							deconvoluted.add(deconvolute(earlyLow, earlyHigh, centerPlusOne, lateLow, lateHigh, tolerance));
						} else if (firstIsOnEarlyBlock&&secondIsOnEarlyBlock) {
							// middle on late block so C
							deconvoluted.add(deconvolute(earlyLow, earlyHigh, center, lateLow, lateHigh, tolerance));
						} else if (!firstIsOnEarlyBlock&&!secondIsOnEarlyBlock) {
							// middle on early block so C+1
							deconvoluted.add(deconvolute(earlyLow, earlyHigh, centerPlusOne, lateLow, lateHigh, tolerance));
						} else {
							// split (both centers are outer edges), so do both
							deconvoluted.add(deconvolute(earlyLow, earlyHigh, center, lateLow, lateHigh, tolerance));
							deconvoluted.add(deconvolute(earlyLow, earlyHigh, centerPlusOne, lateLow, lateHigh, tolerance));
						}
					}
					
					// END THREADABLE SECTION
					
					for (Pair<Stripe, Stripe> pair : deconvoluted) {
						deconvolutedStripes.add(pair.x);
						deconvolutedStripes.add(pair.y);
						
						addRetentionTime(pair.x);
						addRetentionTime(pair.y);
					}
					
					
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
		
		float[] intensities=center.getIntensityArray();
		double[] masses=center.getMassArray();
		
		Range lowerRange=new Range(center.getRange().getStart(), center.getRange().getMiddle());
		Range upperRange=new Range(center.getRange().getMiddle(), center.getRange().getStop());
		
		ArrayList<Peak> lowerPeaks=new ArrayList<Peak>();
		ArrayList<Peak> upperPeaks=new ArrayList<Peak>();
		ArrayList<Peak> removedPeaks=new ArrayList<Peak>();
		
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
			} else {
				removedPeaks.add(new Peak(masses[i], intensities[i]));
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
