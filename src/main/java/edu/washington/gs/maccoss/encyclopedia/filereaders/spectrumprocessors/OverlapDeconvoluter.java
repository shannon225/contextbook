package edu.washington.gs.maccoss.encyclopedia.filereaders.spectrumprocessors;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MSMSBlock;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import gnu.trove.list.array.TFloatArrayList;

public class OverlapDeconvoluter implements SpectrumProcessor {
	private final MassTolerance tolerance;
	private BlockingQueue<MSMSBlock> inputQueue;
	private BlockingQueue<MSMSBlock> outputQueue;
	private final HashMap<Range, TFloatArrayList> retentionTimesByStripe=new HashMap<Range, TFloatArrayList>();
	private final HashMap<Range, TFloatArrayList> ionInjectionTimesByStripe=new HashMap<Range, TFloatArrayList>();
	private final HashMap<Range, TFloatArrayList> truncatedRetentionTimesByStripe=new HashMap<Range, TFloatArrayList>();
	private final HashMap<Range, TFloatArrayList> truncatedIonInjectionTimesByStripe=new HashMap<Range, TFloatArrayList>();

	private Throwable error;

	public OverlapDeconvoluter(MassTolerance tolerance) {
		this.tolerance=tolerance;
	}

	@Override
	public void initialize(StripeFileInterface inputFile, SearchParameters params, BlockingQueue<MSMSBlock> inputQueue, BlockingQueue<MSMSBlock> outputQueue) {
		this.initialize(inputQueue, outputQueue);
	}
	public void initialize(BlockingQueue<MSMSBlock> inputQueue, BlockingQueue<MSMSBlock> outputQueue) {
		this.inputQueue=inputQueue;
		this.outputQueue=outputQueue;
	}

	@Override
	public HashMap<Range, TFloatArrayList> getRetentionTimesByStripe() {
		return retentionTimesByStripe;
	}

	@Override
	public HashMap<Range, TFloatArrayList> getIonInjectionTimesByStripe() {
		return ionInjectionTimesByStripe;
	}

	@Override
	public void run() {
		Range cycleStart=null;
		
		int cycleLength=-1;
		int doubleCycleLength=0;

		// to keep track of whether we're on the early block or the late block
		Range minimumRange=null;
		Range maximumRange=null;

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
						if (cycleStart==null) {
							cycleStart=block.getFragmentScans().get(0).getRange();
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
							doubleCycleLength=cycleLength*2;
							break STARTUP;
						}
						
					}
				}
				
				ArrayList<FragmentScan> deconvolutedStripes=new ArrayList<FragmentScan>();
				BLOCK: for (FragmentScan stripe : block.getFragmentScans()) {
					currentCycle.add(stripe);
					if (currentCycle.size()>doubleCycleLength) {
						currentCycle.removeFirst();
					} else if (currentCycle.size()<doubleCycleLength) {
						continue BLOCK;
					}
					
					FragmentScan earlyLow=null;
					FragmentScan earlyHigh=null;
					FragmentScan center=null;
					FragmentScan lateLow=null;
					FragmentScan lateHigh=null;

					Range target=currentCycle.get(0).getRange();
					float quarterWidth=target.getRange()/4.0f;
					float lowerTarget=target.getMiddle()-quarterWidth;
					float upperTarget=target.getMiddle()+quarterWidth;
					for (int i=1; i<currentCycle.size(); i++) {
						FragmentScan current=currentCycle.get(i);
						Range range=current.getRange();
						if (range.equals(target)) {
							center=current;
						} else if (range.contains(lowerTarget)) {
							if (center==null) {
								earlyLow=current;
							} else if (lateLow==null) {
								lateLow=current;
							}
						} else if (range.contains(upperTarget)) {
							if (center==null) {
								earlyHigh=current;
							} else if (lateHigh==null) {
								lateHigh=current;
							}
						}
						if (lateLow!=null&&lateHigh!=null) {
							break;
						}
					}

					try {
						Pair<FragmentScan, FragmentScan> pair=deconvolute(earlyLow, earlyHigh, center, lateLow, lateHigh, tolerance);
						deconvolutedStripes.add(pair.x);
						deconvolutedStripes.add(pair.y);
						
						addRetentionTime(pair.x);
						addRetentionTime(pair.y);
					} catch (Exception e) {
						Logger.errorLine("Error deconvoluting sample!");
						Logger.errorException(e);
					}

				}
				outputQueue.put(new MSMSBlock(block.getPrecursorScans(), deconvolutedStripes));
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("DIA writing interrupted!");
			Logger.errorException(ie);
		} catch (Throwable t) {
			Logger.errorLine("Overlap deconvolution failed!");
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
	
	protected static Pair<FragmentScan, FragmentScan> deconvolute(FragmentScan earlyLow, FragmentScan earlyHigh, FragmentScan center, FragmentScan lateLow, FragmentScan lateHigh, MassTolerance tolerance) {
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
		
		FragmentScan lowerStripe=getDeconvolutedStripe(center, lowerRange, lowerPeaks, false);
		FragmentScan upperStripe=getDeconvolutedStripe(center, upperRange, upperPeaks, true);
		
		Pair<FragmentScan, FragmentScan> deconvoluted=new Pair<FragmentScan, FragmentScan>(lowerStripe, upperStripe);
		return deconvoluted;
	}

	private static FragmentScan getDeconvolutedStripe(FragmentScan center, Range lowerRange, ArrayList<Peak> lowerPeaks, boolean useNegativeScanNumber) {
		Pair<double[], float[]> arrays=Peak.toArrays(lowerPeaks);
		int scanNumber=useNegativeScanNumber?(Integer.MAX_VALUE-center.getSpectrumIndex()):center.getSpectrumIndex();
		
		FragmentScan lowerStripe=new FragmentScan(center.getSpectrumName(), center.getPrecursorName(), scanNumber, center.getScanStartTime(), center.getFraction(), center.getIonInjectionTime(), lowerRange.getStart(), lowerRange.getStop(), arrays.x, arrays.y);
		return lowerStripe;
	}

	private static float getIntensity(MassTolerance tolerance, double mass, FragmentScan stripe) {
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
