package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PecanLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import gnu.trove.map.hash.TDoubleFloatHashMap;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.procedure.TDoubleFloatProcedure;

public class PecanOneFragmentationModel extends AbstractPecanFragmentationModel {
	private final String accession;
	public PecanOneFragmentationModel(FastaEntry peptide, AminoAcidConstants aaConstants) {
		super(peptide.getSequence(), aaConstants);
		this.accession=peptide.getAccession();
	}
	public PecanLibraryEntry getUnitSpectrum(byte precursorCharge, SearchParameters params) {
		double[] ions=getPrimaryIons(params.getFragType(), precursorCharge);
		float[] unitIntensities=new float[ions.length];
		Arrays.fill(unitIntensities, 1.0f);
		float euclidianDistance=(float)Math.sqrt(ions.length);

		String sequence=getModifiedSequence();
		double precursorMZ=params.getAAConstants().getChargedMass(sequence, precursorCharge);

		return new PecanLibraryEntry(accession, precursorMZ, precursorCharge, sequence, 1, 0.0f, 0, ions, unitIntensities, false, euclidianDistance);
	}
	public PecanLibraryEntry getPecanSpectrum(byte precursorCharge, double[] sortedBinCounterKeys, TDoubleIntHashMap binCounter, Range fragmentationRange, SearchParameters params, boolean isDecoy) {
		TDoubleFloatHashMap peakMap=new TDoubleFloatHashMap();
		double[] ions=getPrimaryIons(params.getFragType(), precursorCharge);
		float totalOfSquares=0.0f;
		for (int i=0; i<ions.length; i++) {
			if (!fragmentationRange.contains((float)ions[i])) {
				continue;
			}
			double[] matches=params.getFragmentTolerance().getMatches(sortedBinCounterKeys, ions[i]);
			
			int total=0;
			if (matches.length>0) {
				for (int j=0; j<matches.length; j++) {
					total+=binCounter.get(matches[j]);
				}
			}
			if (total==0) total=1; // TODO not quite adding one pseudocount (which would probably be more robust)
			
			float score=100.0f/total;
			peakMap.put(ions[i], score);
			totalOfSquares+=score*score;
		}
		
		final float euclidianDistance=(float)Math.sqrt(totalOfSquares);
		final ArrayList<Peak> peaks=new ArrayList<Peak>();
		
		peakMap.forEachEntry(new TDoubleFloatProcedure() {
			public boolean execute(double arg0, float arg1) {
				peaks.add(new Peak(arg0, arg1/euclidianDistance));
				return true;
			}
		});
		
		Collections.sort(peaks);
		Pair<double[], float[]> arrays=Peak.toArrays(peaks);
		
		String sequence=getModifiedSequence();
		double precursorMZ=params.getAAConstants().getChargedMass(sequence, precursorCharge);

		return new PecanLibraryEntry(accession, precursorMZ, precursorCharge, sequence, 1, 0.0f, 0, arrays.x, arrays.y, isDecoy, euclidianDistance);	
	}
}
