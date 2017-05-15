package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.procedure.TDoubleIntProcedure;
import gnu.trove.set.hash.TDoubleHashSet;

/**
 * follows BackgroundGenerator from Pecan, but deviates in how it is generated/used
 * @author searleb
 *
 */
public class BackgroundFrequencyCalculator implements BackgroundFrequencyInterface {
	private static final int MASS_COUNTER_BIN_COUNT=1000;
	private final double[] binBoundaries;
	private final TDoubleIntHashMap[] binCounters;
	private final double[][] sortedMapKeys;
	private final int[] numberOfSpectra;
	
	BackgroundFrequencyCalculator(double[] binBoundaries, TDoubleIntHashMap[] binCounters, int[] numberOfSpectra) {
		this.binBoundaries=binBoundaries;
		this.binCounters=binCounters;
		sortedMapKeys=new double[binCounters.length][];
		for (int i=0; i<binCounters.length; i++) {
			sortedMapKeys[i]=binCounters[i].keys();
			Arrays.sort(sortedMapKeys[i]);
		}
		this.numberOfSpectra=numberOfSpectra;
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyInterface#getRoundedMassCounters(double, edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance)
	 */
	@Override
	public Pair<double[], float[]> getRoundedMassCounters(double precursorMz, MassTolerance tolerance) {
		final double[] bestMass=new double[MASS_COUNTER_BIN_COUNT];
		int binIndex=getBinIndex(precursorMz);
		if (binIndex<0||binIndex>=binCounters.length) {
			return new Pair<double[], float[]>(bestMass, new float[MASS_COUNTER_BIN_COUNT]);
		}
		TDoubleIntHashMap counter=binCounters[binIndex];

		final int[] bestMassCount=new int[MASS_COUNTER_BIN_COUNT];
		Arrays.fill(bestMassCount, Integer.MIN_VALUE);
		counter.forEachEntry(new TDoubleIntProcedure() {
			@Override
			public boolean execute(double a, int b) {
				int index=(int)a; // truncates! Should be fine up to 1000 m/z
				if (index<MASS_COUNTER_BIN_COUNT&&b>bestMassCount[index]) {
					bestMass[index]=a;
					bestMassCount[index]=b;
				}
				return true;
			}
		});
		return new Pair<double[], float[]>(bestMass, getFrequencies(bestMass, precursorMz, tolerance));
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
		
		TDoubleIntHashMap[] binCounters=new TDoubleIntHashMap[binBoundaries.length-1];
		int[] numberOfSpectra=new int[binBoundaries.length-1];
		for (int i=0; i<binCounters.length; i++) {
			binCounters[i]=new TDoubleIntHashMap();
			numberOfSpectra[i]=1; // add pseudocount
		}

		Logger.logLine("Creating background frequency calculator from "+diafile.getOriginalFileName()+"...");
		for (Range range : ranges) {
			double targetMz=range.getMiddle();
			Logger.logLine("Processing "+range.toString()+" m/z");
			int index=Arrays.binarySearch(binBoundaries, targetMz);
			if (index<0) {
				index=(-(index+1))-1;
			}
			if (index<0||index>=binCounters.length) continue;
			
			ArrayList<Stripe> stripes=diafile.getStripes(targetMz, -Float.MAX_VALUE, Float.MAX_VALUE, false);
			for (Stripe stripe : stripes) {
				double[] ions=stripe.getMassArray();
				
				numberOfSpectra[index]++;
				for (double ion : ions) {
					binCounters[index].adjustOrPutValue(ion, 1, 1);
				}
			}
		}

		return new BackgroundFrequencyCalculator(binBoundaries, binCounters, numberOfSpectra);
	}
	
	public static BackgroundFrequencyInterface generateBackground(StripeFileInterface diafile, LibraryInterface library) throws DataFormatException, SQLException, IOException {
		TDoubleHashSet boundaries=new TDoubleHashSet();
		for (Range range : diafile.getRanges().keySet()) {
			boundaries.add(range.getStart());
			boundaries.add(range.getStop());
		}
		double[] binBoundaries=boundaries.toArray();
		Arrays.sort(binBoundaries);
		
		TDoubleIntHashMap[] binCounters=new TDoubleIntHashMap[binBoundaries.length-1];
		int[] numberOfSpectra=new int[binBoundaries.length-1];
		for (int i=0; i<binCounters.length; i++) {
			binCounters[i]=new TDoubleIntHashMap();
			numberOfSpectra[i]=1; // add pseudocount
		}

		if (library!=null) {
			ArrayList<LibraryEntry> allEntries=library.getAllEntries(false);
			for (LibraryEntry entry : allEntries) {
				double[] ions=entry.getMassArray();
				int index=Arrays.binarySearch(binBoundaries, entry.getPrecursorMZ());
				index=(-(index+1))-1;

				if (index<0||index>=binCounters.length) continue;

				numberOfSpectra[index]++;
				for (double ion : ions) {
					binCounters[index].adjustOrPutValue(ion, 1, 1);
				}
			}
		}

		return new BackgroundFrequencyCalculator(binBoundaries, binCounters, numberOfSpectra);
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyInterface#getFrequencies(double[], double, edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance)
	 */
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
			double[] matches=tolerance.getMatches(sortedMapKeys[index], ions[i]);
			
			if (matches.length>0) {
				for (int j=0; j<matches.length; j++) {
					counters[i]+=binCounters[index].get(matches[j]);
				}
			}
		}
		return getFrequencies(counters, numberOfLibraryEntries);
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
	
	private int getBinIndex(double precursorMz) {
		int index=Arrays.binarySearch(binBoundaries, precursorMz);
		if (index<0) {
			index=(-(index+1))-1;
		}
		return index;
	}
	
	private float[] getFrequencies(int[] counters, int norm) {
		float[] frequencies=new float[counters.length];
		for (int i=0; i<frequencies.length; i++) {
			frequencies[i]=counters[i]/(float)norm;//numberOfLibraryEntries;
		}
		return frequencies;
	}
}
