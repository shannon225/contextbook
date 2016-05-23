package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.IOException;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.map.hash.TDoubleIntHashMap;
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
	private final int totalCount;
	
	BackgroundFrequencyCalculator(double[] binBoundaries, TDoubleIntHashMap[] binCounters) {
		this.binBoundaries=binBoundaries;
		this.binCounters=binCounters;
		sortedMapKeys=new double[binCounters.length][];
		int total=0;
		for (int i=0; i<binCounters.length; i++) {
			total+=General.sum(binCounters[i].values());
			sortedMapKeys[i]=binCounters[i].keys();
			Arrays.sort(sortedMapKeys[i]);
		}
		this.totalCount=total;
	}
	
	public static BackgroundFrequencyCalculator generateBackground(StripeFile diafile, LibraryFile library) throws DataFormatException, SQLException, IOException {
		TDoubleHashSet boundaries=new TDoubleHashSet();
		for (Range range : diafile.getRanges().keySet()) {
			boundaries.add(range.getStart());
			boundaries.add(range.getStop());
		}
		double[] binBoundaries=boundaries.toArray();
		Arrays.sort(binBoundaries);
		
		TDoubleIntHashMap[] binCounters=new TDoubleIntHashMap[binBoundaries.length-1];
		for (int i=0; i<binCounters.length; i++) {
			binCounters[i]=new TDoubleIntHashMap();
		}
		
		for (LibraryEntry entry : library.getAllEntries(false)) {
			double[] ions=entry.getMassArray();
			int index=Arrays.binarySearch(binBoundaries, entry.getPrecursorMZ());
			index=(-(index+1))-1;

			if (index<0||index>=binCounters.length) continue;
			
			for (double ion : ions) {
				binCounters[index].adjustOrPutValue(ion, 1, 1);
			}
		}

		return new BackgroundFrequencyCalculator(binBoundaries, binCounters);
	}
	
	public float[] getFrequencies(double[] ions, double precursorMz, MassTolerance tolerance) {
		int[] counters=new int[ions.length];
		Arrays.fill(counters, 1); // add pseudocount
		
		int index=Arrays.binarySearch(binBoundaries, precursorMz);
		index=(-(index+1))-1;
		if (index<0||index>=binCounters.length) {
			return getFrequencies(counters);
		}
		
		for (int i=0; i<ions.length; i++) {
			double[] matches=tolerance.getMatches(sortedMapKeys[index], ions[i]);
			
			if (matches.length>0) {
				for (int j=0; j<matches.length; j++) {
					counters[i]+=binCounters[index].get(matches[j]);
				}
			}
		}
		return getFrequencies(counters);
	}
	
	private float[] getFrequencies(int[] counters) {
		float[] frequencies=new float[counters.length];
		for (int i=0; i<frequencies.length; i++) {
			frequencies[i]=counters[i]/totalCount;
		}
		return frequencies;
	}
}
