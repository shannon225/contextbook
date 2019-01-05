package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.set.hash.TDoubleHashSet;

public class BinnedBackgroundFrequencyCalculator implements BackgroundFrequencyInterface {
	private static final int FRACTIONAL_INCREMENT=2000; // 5 ppm at 100
	private static final int LARGEST_MASS=1600;
	private static final int BIN_LENGTH=LARGEST_MASS*FRACTIONAL_INCREMENT;
	private final double[] binBoundaries;
	private final int[][] binCounters;
	private final int[] numberOfSpectra;

	public BinnedBackgroundFrequencyCalculator(double[] binBoundaries, int[][] binCounters, int[] numberOfSpectra) {
		this.binBoundaries=binBoundaries;
		this.binCounters=binCounters;
		this.numberOfSpectra=numberOfSpectra;
	}

	@Override
	public LibraryBackgroundInterface getLibraryBackground(final double precursorMz, final MassTolerance tolerance) {
		return new LibraryBackgroundInterface() {
			@Override
			public float getFraction(double mass) {
				return getFrequencies(new double[] {mass}, precursorMz, tolerance)[0];
			}
		};
	}
	
	@Override
	public float[] getFrequencies(double[] ions, double precursorMz, MassTolerance tolerance) {
		int[] counters=new int[ions.length];
		Arrays.fill(counters, 1); // add pseudocount
		
		int index=getBinIndex(precursorMz);
		if (index<0||index>=binCounters.length) {
			return getFrequencies(counters, 1);
		}
		
		int numberOfLibraryEntries=numberOfSpectra[index];
		
		for (int i=0; i<ions.length; i++) {
			int binIndex=getIndex(ions[i]);
			counters[i]=binCounters[index][binIndex];

			int lowerIndex=binIndex-1;
			while (tolerance.equals(ions[i], getMass(lowerIndex))) {
				counters[i]+=binCounters[index][lowerIndex];
				lowerIndex--;
			}
			int upperIndex=binIndex+1;
			while (tolerance.equals(ions[i], getMass(upperIndex))) {
				counters[i]+=binCounters[index][upperIndex];
				upperIndex++;
			}
		}
		return getFrequencies(counters, numberOfLibraryEntries);
	}

	@Override
	public Pair<double[], float[]> getRoundedMassCounters(double precursorMz, MassTolerance tolerance) {
		final double[] bestMass=new double[LARGEST_MASS];
		int binIndex=getBinIndex(precursorMz);
		if (binIndex<0||binIndex>=binCounters.length) {
			return new Pair<double[], float[]>(bestMass, new float[LARGEST_MASS]);
		}

		int[] bestFrequencyCounter=new int[LARGEST_MASS];
		for (int i=0; i<binCounters[binIndex].length; i++) {
			float mass=i/(float)FRACTIONAL_INCREMENT;
			int index=(int)mass;
			if (bestFrequencyCounter[index]<binCounters[binIndex][i]) {
				bestFrequencyCounter[index]=binCounters[binIndex][i];
				bestMass[index]=mass;
			}
		}

		return new Pair<double[], float[]>(bestMass, getFrequencies(bestMass, precursorMz, tolerance));
	}
	
	int getNumberOfLibraryEntries(double precursorMz) {
		int index=getBinIndex(precursorMz);
		int numberOfLibraryEntries=numberOfSpectra[index];
		return numberOfLibraryEntries;
	}

	private int getBinIndex(double precursorMz) {
		int index=Arrays.binarySearch(binBoundaries, precursorMz);
		if (index<0) {
			index=(-(index+1))-1;
		}
		return index;
	}

	public static BackgroundFrequencyInterface generateBackground(StripeFileInterface diafile) throws DataFormatException, SQLException, IOException {
		TDoubleHashSet boundaries=new TDoubleHashSet();
		ArrayList<Range> ranges=new ArrayList<Range>(diafile.getRanges().keySet());
		Collections.sort(ranges);
		
		for (Range range : ranges) {
			boundaries.add(range.getStart());
			boundaries.add(range.getStop());
		}
		double[] binBoundaries=boundaries.toArray();
		Arrays.sort(binBoundaries);
		
		int[][] binCounters=new int[binBoundaries.length-1][];
		int[] numberOfSpectra=new int[binBoundaries.length-1];
		for (int i=0; i<binCounters.length; i++) {
			binCounters[i]=new int[BIN_LENGTH];
			numberOfSpectra[i]=1; // add pseudocount
		}

		for (Range range : ranges) {
			double targetMz=range.getMiddle();
			Logger.logLine("Processing "+range.toString()+" m/z");
			int index=Arrays.binarySearch(binBoundaries, targetMz);
			if (index<0) {
				index=(-(index+1))-1;
			}
			if (index<0||index>=binCounters.length) continue;
			
			ArrayList<FragmentScan> stripes=diafile.getStripes(targetMz, -Float.MAX_VALUE, Float.MAX_VALUE, false);
			for (FragmentScan stripe : stripes) {
				double[] ions=stripe.getMassArray();

				numberOfSpectra[index]++;
				for (double ion : ions) {
					binCounters[index][getIndex(ion)]++;
				}
			}
		}

		return new BinnedBackgroundFrequencyCalculator(binBoundaries, binCounters, numberOfSpectra);
	}
	
	static int getIndex(double mass) {
		int index=(int)(mass*FRACTIONAL_INCREMENT); // truncates
		if (index<0) return 0;
		if (index>=BIN_LENGTH) return BIN_LENGTH-1;
		return index;
	}

	static double getMass(int binIndex) {
		return binIndex/(double)FRACTIONAL_INCREMENT;
	}
	
	private static float[] getFrequencies(int[] counters, int norm) {
		float[] frequencies=new float[counters.length];
		for (int i=0; i<frequencies.length; i++) {
			frequencies[i]=counters[i]/(float)norm;
		}
		return frequencies;
	}
}
