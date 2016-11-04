package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.procedure.TDoubleIntProcedure;
import gnu.trove.set.hash.TDoubleHashSet;

/**
 * follows BackgroundGenerator from Pecan, but deviates in how it is generated/used
 * @author searleb
 *
 */
public class BackgroundFrequencyCalculator {
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
	
	/**
	 * for display only
	 * @param precursorMz
	 * @return
	 */
	public int[] getRoundedMassCounters() {
		final TDoubleIntHashMap giantCounter=new TDoubleIntHashMap();
		for (TDoubleIntHashMap binCounter : binCounters) {
			binCounter.forEachEntry(new TDoubleIntProcedure() {
				@Override
				public boolean execute(double a, int b) {
					giantCounter.adjustOrPutValue(a, b, b);
					return true;
				}
			});
		}

		final double[] bestMass=new double[1000];
		final int[] bestMassCount=new int[1000];
		giantCounter.forEachEntry(new TDoubleIntProcedure() {
			@Override
			public boolean execute(double a, int b) {
				int index=(int)a; // truncates! Should be fine up to 1000 m/z
				if (b>bestMassCount[index]) {
					bestMass[index]=a;
					bestMassCount[index]=b;
				}
				return false;
			}
		});

		int[] totalCounters=new int[1000];
		
		for (int index=0; index<sortedMapKeys.length; index++) {
			int[] counters=new int[1000];
			double[] matches=sortedMapKeys[index];
			for (int j=0; j<matches.length; j++) {
				int count=binCounters[index].get(matches[j]);
				int i=(int) Math.round(matches[j]);
				if (i<counters.length) {
					if (counters[i]<count) {
						counters[i]=count;
					}
				}
			}
			for (int i=0; i<totalCounters.length; i++) {
				totalCounters[i]+=counters[i];
			}
		}
		return totalCounters;
	}
	
	public static BackgroundFrequencyCalculator generateBackground(StripeFileInterface diafile, LibraryInterface library) throws DataFormatException, SQLException, IOException {
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
	
	public float[] getFrequencies(double[] ions, double precursorMz, MassTolerance tolerance) {
		int[] counters=new int[ions.length];
		Arrays.fill(counters, 1); // add pseudocount
		
		int index=Arrays.binarySearch(binBoundaries, precursorMz);
		index=(-(index+1))-1;
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
	
	private float[] getFrequencies(int[] counters, int norm) {
		float[] frequencies=new float[counters.length];
		for (int i=0; i<frequencies.length; i++) {
			frequencies[i]=counters[i]/(float)norm;//numberOfLibraryEntries;
		}
		return frequencies;
	}
}
