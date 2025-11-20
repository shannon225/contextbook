package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;
import java.util.Optional;
import java.util.TreeMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;

public class SpectrumPeakFilter {
	protected final static PeakIntensityComparator intensityComparator=PeakIntensityComparator.DEFAULT_INTENSITY_COMPARATOR;
	protected final static int NUM_PEAKS_PER_BIN=10;
	protected final static double BIN_SIZE=20.0;// m/z
	protected final static int NUM_OF_BINS=100; // only consider up to 2,000 m/z
												// (anything over will be placed
												// in the last bin)

	public static FragmentScan filterPeaks(FragmentScan stripe) {
		return filterPeaks(stripe, BIN_SIZE, NUM_OF_BINS, NUM_PEAKS_PER_BIN);
	}
	

	public static PrecursorScan combineAndFilterPeaks(PrecursorScan precursor, MassTolerance tolerance) {
		ArrayList<Peak> peaks=Peak.fromArrays(precursor.getMassArray(), precursor.getIntensityArray(), precursor.getIonMobilityArray());
		
		peaks=combinePeaks(peaks, tolerance);
		peaks=filterPeaks(peaks, BIN_SIZE, NUM_OF_BINS, NUM_PEAKS_PER_BIN);
		Triplet<double[], float[], Optional<float[]>> arrays=Peak.toArrays(peaks);
		return new PrecursorScan(precursor.getSpectrumName(), precursor.getSpectrumIndex(), precursor.getScanStartTime(), precursor.getFraction(), precursor.getIsolationWindowLower(), precursor.getIsolationWindowUpper(), precursor.getIonInjectionTime(), 
				arrays.x, arrays.y, arrays.z);
	}

	public static FragmentScan combineAndFilterPeaks(FragmentScan stripe, MassTolerance tolerance) {
		ArrayList<Peak> peaks=Peak.fromArrays(stripe.getMassArray(), stripe.getIntensityArray(), stripe.getIonMobilityArray());
		
		peaks=combinePeaks(peaks, tolerance);
		peaks=filterPeaks(peaks, BIN_SIZE, NUM_OF_BINS, NUM_PEAKS_PER_BIN);
		Triplet<double[], float[], Optional<float[]>> arrays=Peak.toArrays(peaks);

		return new FragmentScan(stripe.getSpectrumName(), stripe.getPrecursorName(), stripe.getSpectrumIndex(), stripe.getScanStartTime(), stripe.getFraction(), stripe.getIonInjectionTime(), stripe.getIsolationWindowLower(), stripe.getIsolationWindowUpper(),
				arrays.x, arrays.y, arrays.z, stripe.getPrecursorCharge());
	}

	public static FragmentScan filterPeaks(FragmentScan stripe, double binSize, int numBins, int numPeaksPerBin) {
		ArrayList<Peak> peaks=Peak.fromArrays(stripe.getMassArray(), stripe.getIntensityArray(), stripe.getIonMobilityArray());
		peaks=filterPeaks(peaks, binSize, numBins, numPeaksPerBin);
		Triplet<double[], float[], Optional<float[]>> arrays=Peak.toArrays(peaks);

		return new FragmentScan(stripe.getSpectrumName(), stripe.getPrecursorName(), stripe.getSpectrumIndex(), stripe.getScanStartTime(), stripe.getFraction(), stripe.getIonInjectionTime(), stripe.getIsolationWindowLower(), stripe.getIsolationWindowUpper(),
				arrays.x, arrays.y, arrays.z, stripe.getPrecursorCharge());
	}

	public static ArrayList<Peak> filterPeaks(ArrayList<Peak> peaks, double binSize, int numBins, int numPeaksPerBin) {
		@SuppressWarnings("unchecked")
		
		ArrayList<Peak>[] bins=new ArrayList[numBins];
		for (int i=0; i<bins.length; i++) {
			bins[i]=new ArrayList<Peak>();
		}

		for (Peak peak : peaks) {
			int index=getIndex(peak.mass, binSize, numBins);
			bins[index].add(peak);
		}

		ArrayList<Peak> filtered=new ArrayList<>();
		for (ArrayList<Peak> list : bins) {
			Collections.sort(list, intensityComparator);

			int stopIndex=list.size()-numPeaksPerBin;
			if (stopIndex>0) {
				for (int i=list.size()-1; i>=stopIndex; i--) {
					filtered.add(list.get(i));
				}
			} else {
				// just add all peaks
				filtered.addAll(list);
			}
		}

		// final sort on mass
		Collections.sort(filtered);
		return filtered;
	}
	
	private static int getIndex(double mz, double binSize, int numBins) {
		int index=(int)(mz/binSize);
		if (index>=numBins) return numBins-1;
		if (index<0) return 0;
		return index;
	}
	
	public static ArrayList<Peak> combinePeaks(ArrayList<Peak> peaks, final MassTolerance tolerance) {
		TreeMap<Double, Peak> massSortedPeaks=new TreeMap<Double, Peak>();

		ArrayList<Peak> intensitySortedPeaks=new ArrayList<Peak>(peaks);
		Collections.sort(intensitySortedPeaks, intensityComparator);
		
		for (Peak peak : intensitySortedPeaks) {
			Map.Entry<Double, Peak> lo=massSortedPeaks.floorEntry(peak.getMass());
			Map.Entry<Double, Peak> hi=massSortedPeaks.ceilingEntry(peak.getMass());

			Peak best=null;
			if (lo!=null&&tolerance.equals(lo.getValue().getMass(), peak.getMass())) best=lo.getValue();
			if (hi!=null&&tolerance.equals(hi.getValue().getMass(), peak.getMass())) {
				if (best==null||hi.getValue().getIntensity()>best.getIntensity()) best=hi.getValue();
			}

			if (best==null) {
				massSortedPeaks.put(peak.getMass(), peak);
			} else {
				Peak merged=Peak.mergePeaks(best, peak);
				massSortedPeaks.remove(best.getMass());
				massSortedPeaks.put(merged.getMass(), merged);
			}
		}
		return new ArrayList<Peak>(massSortedPeaks.values());
	}
}
