package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
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
		int index=Arrays.binarySearch(rts, rt);
		if (index<0) {
			// insertion point
			index=-(index+1);
			
			// if we're not at the first bin, choose the previous (precursor) scan
			if (index>0) index--;
			if (index>=rts.length) index--;
		}
		
		PrecursorScan scan=precursors.get(index);
		float[] intensities=scan.getIntensityArray();
		double[] masses=scan.getMassArray();
		
		Peak[] isotopeIntensities=new Peak[isotopes.length];
		for (int i=0; i<isotopes.length; i++) {
			byte isotope=isotopes[i];
			
			double target=mz+isotope*charge*MassConstants.neutronMass;
			int[] indicies=tolerance.getIndicies(masses, target);
			float intensity=0.0f;
			double weightedMz=0.0;
			for (int j=0; j<indicies.length; j++) {
				intensity+=intensities[indicies[j]];
				weightedMz+=intensities[indicies[j]]*masses[indicies[j]];
			}
			isotopeIntensities[i]=new Peak(weightedMz/intensity, intensity);
		}
		return isotopeIntensities;
	}
}
