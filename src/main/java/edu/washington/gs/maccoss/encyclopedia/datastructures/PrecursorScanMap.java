package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import gnu.trove.list.array.TFloatArrayList;

//@Immutable
public class PrecursorScanMap {
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
	
	public Peak[] getIsotopePacket(double mz, float rt, byte charge, MassTolerance tolerance) {
		if (rts.length==0) return new Peak[0];
		
		PrecursorScan scan=getNearestScan(mz, rt);
		if (scan==null) return new Peak[0];
		
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
