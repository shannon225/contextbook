package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;

public class SpectrumPeakFilter {
	private final static PeakIntensityComparator intensityComparator=new PeakIntensityComparator();
	private final static int NUM_PEAKS_PER_BIN=10;
	private final static double BIN_SIZE=20.0;// m/z
	private final static int NUM_OF_BINS=100; // only consider up to 2,000 m/z
												// (anything over will be placed
												// in the last bin)

	public static FragmentScan filterPeaks(FragmentScan stripe) {
		ArrayList<PeakChromatogram> peaks=PeakChromatogram.fromChromatogramArrays(stripe.getMassArray(), stripe.getIntensityArray(), new float[stripe.getMassArray().length]);
		peaks=filterPeaks(peaks);
		Triplet<double[], float[], float[]> arrays=PeakChromatogram.toChromatogramArrays(peaks);

		return new FragmentScan(stripe.getSpectrumName(), stripe.getPrecursorName(), stripe.getSpectrumIndex(), stripe.getScanStartTime(), stripe.getFraction(), stripe.getIonInjectionTime(), stripe.getIsolationWindowLower(), stripe.getIsolationWindowUpper(),
				arrays.x, arrays.y, stripe.getCharge());
	}

	public static LibraryEntry filterPeaks(LibraryEntry entry, AminoAcidConstants aaConstants) {
		ArrayList<PeakChromatogram> peaks=PeakChromatogram.fromChromatogramArrays(entry.getMassArray(), entry.getIntensityArray(), entry.getCorrelationArray());
		peaks=filterPeaks(peaks);
		Triplet<double[], float[], float[]> arrays=PeakChromatogram.toChromatogramArrays(peaks);

		return new LibraryEntry(entry.getSource(), entry.getAccessions(), entry.getSpectrumIndex(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), entry.getPeptideModSeq(), entry.getCopies(),
				entry.getRetentionTime(), entry.getScore(), arrays.x, arrays.y, arrays.z, aaConstants);
	}

	public static ArrayList<PeakChromatogram> filterPeaks(ArrayList<PeakChromatogram> peaks) {
		@SuppressWarnings("unchecked")
		ArrayList<PeakChromatogram>[] bins=new ArrayList[NUM_OF_BINS];
		for (int i=0; i<bins.length; i++) {
			bins[i]=new ArrayList<PeakChromatogram>();
		}

		for (PeakChromatogram peak : peaks) {
			int index=getIndex(peak.mass);
			bins[index].add(peak);
		}

		ArrayList<PeakChromatogram> filtered=new ArrayList<>();
		for (ArrayList<PeakChromatogram> list : bins) {
			Collections.sort(list, intensityComparator);

			int stopIndex=list.size()-NUM_PEAKS_PER_BIN;
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

	private static int getIndex(double mz) {
		int index=(int)(mz/BIN_SIZE);
		if (index>=NUM_OF_BINS) return NUM_OF_BINS-1;
		if (index<0) return 0;
		return index;
	}
}
