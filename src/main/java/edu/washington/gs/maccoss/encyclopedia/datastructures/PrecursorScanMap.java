package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;

//@Immutable
public class PrecursorScanMap {
	public static final int monoisotopicIndex=1;
	public static final int plusOneIndex=2;
	public static final byte[] isotopes=new byte[] {-1, 0, 1, 2};
	
	private final float[] rts;
	private final ArrayList<PrecursorScan> precursors;

	public PrecursorScanMap(ArrayList<PrecursorScan> precursors) {
		this.precursors=precursors;
		Collections.sort(this.precursors);
		TFloatArrayList rts=new TFloatArrayList();
		for (PrecursorScan scan : precursors) {
			rts.add(scan.getScanStartTime());
		}
		this.rts=rts.toArray();
	}
	
	public Pair<float[], TFloatArrayList[]> integrateIsotopePacket(double mz, Range rtRange, byte charge, MassTolerance tolerance) {
		float[] summedIntensity=new float[isotopes.length];
		TFloatArrayList[] deltaMasses=new TFloatArrayList[isotopes.length];
		for (int i = 0; i < deltaMasses.length; i++) {
			deltaMasses[i]=new TFloatArrayList();
		}
		
		int index=Arrays.binarySearch(rts, rtRange.getStart());
		if (index<0) {
			// insertion point
			index=-(index+1);
		}
		if (index<0) index=0;
		if (index>=precursors.size()) return new Pair<float[], TFloatArrayList[]>(summedIntensity, deltaMasses);
		
		Peak[] previousIsotopeIntensities = null;
		
		// integrate exclusively (only scans within the range
		for (int i = index+1; i < rts.length; i++) {
			if (rtRange.contains(rts[i])) {
				Peak[] isotopeIntensities = getPacket(precursors.get(i), mz, charge, tolerance);
				if (previousIsotopeIntensities==null) {
					previousIsotopeIntensities=isotopeIntensities;
					continue;
				}
				
				float deltaTime=rts[i]-rts[i-1];
				
				for (int j = 0; j < isotopeIntensities.length; j++) {
					float area=deltaTime*(isotopeIntensities[j].getIntensity()+previousIsotopeIntensities[j].getIntensity())/2.0f;
					summedIntensity[j]+=area;
					
					double target=mz+(isotopes[j]*MassConstants.neutronMass/charge);
					double deltaMass=isotopeIntensities[j].getMass()-target;
					if (tolerance.isRelativeTolerance()) {
						deltaMasses[j].add((float)(1000000.0*deltaMass/target));
					} else {
						deltaMasses[j].add((float)deltaMass);
					}
				}
				isotopeIntensities=previousIsotopeIntensities;
			} else {
				// outside range
				break;
			}
		}
		
		if (General.sum(summedIntensity)==0.0f&&previousIsotopeIntensities!=null) {
			// in case there's only one point
			for (int j = 0; j < previousIsotopeIntensities.length; j++) {
				summedIntensity[j]=previousIsotopeIntensities[j].intensity;
			}
		}
		
		return new Pair<float[], TFloatArrayList[]>(summedIntensity, deltaMasses);
	}
	
	public Peak[] getIsotopePacket(double mz, float rt, byte charge, MassTolerance tolerance) {
		if (rts.length==0) return new Peak[0];
		
		PrecursorScan scan=getNearestScan(mz, rt);
		if (scan==null) return new Peak[0];
		
		Peak[] isotopeIntensities = getPacket(scan, mz, charge, tolerance);
		return isotopeIntensities;
	}

	private Peak[] getPacket(PrecursorScan scan, double mz, byte charge, MassTolerance tolerance) {
		float[] intensities=scan.getIntensityArray();
		double[] masses=scan.getMassArray();
		
		Peak[] isotopeIntensities=new Peak[isotopes.length];
		for (int i=0; i<isotopes.length; i++) {
			byte isotope=isotopes[i];
			
			double target=mz+(isotope*MassConstants.neutronMass/charge);
			int[] indicies=tolerance.getIndicies(masses, target);
			float intensity=0.0f;
			double bestMz=target;
			float bestPeakIntensity=0.0f;
			for (int j=0; j<indicies.length; j++) {
				intensity+=intensities[indicies[j]];
				if (intensities[indicies[j]]>bestPeakIntensity) {
					bestPeakIntensity=intensities[indicies[j]];
					bestMz=masses[indicies[j]];
				}
				
			}
			
			isotopeIntensities[i]=new Peak(bestMz, intensity);
		}
		return isotopeIntensities;
	}
	
	private PrecursorScan getNearestScan(double mz, float rt) {
		int index=Arrays.binarySearch(rts, rt);
		if (index<0) {
			// insertion point
			index=-(index+1);
		}
		
		int upperIndex=index;
		PrecursorScan scanAbove=null;
		while (upperIndex<rts.length) {
			scanAbove=precursors.get(upperIndex);

			if (mz>scanAbove.getIsolationWindowLower()&&mz<scanAbove.getIsolationWindowUpper()) {
				break;
			} else {
				upperIndex++;
			}
		}
		float deltaRTAbove=scanAbove==null?Float.MAX_VALUE:Math.abs(rt-scanAbove.getScanStartTime());
		
		int lowerIndex=index-1;
		PrecursorScan scanBelow=null;
		while (lowerIndex>=0) {
			scanBelow=precursors.get(lowerIndex);

			if (mz>scanBelow.getIsolationWindowLower()&&mz<scanBelow.getIsolationWindowUpper()) {
				break;
			} else {
				lowerIndex--;
			}
		}
		float deltaRTBelow=scanBelow==null?Float.MAX_VALUE:Math.abs(rt-scanBelow.getScanStartTime());
		
		if (deltaRTBelow>deltaRTAbove) {
			return scanAbove;
		} else {
			return scanBelow;
		}
	}
}
