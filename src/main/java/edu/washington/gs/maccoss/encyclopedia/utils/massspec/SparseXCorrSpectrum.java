package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.utils.SparseIndexMap;
import gnu.trove.procedure.TIntObjectProcedure;

public class SparseXCorrSpectrum {
	private final float fragmentBinSize;
	private final int[] indices;
	private final double[] masses;
	private final float[] intensities;
	private final int length;
	private final double precursorMz;
	
	SparseXCorrSpectrum(SparseIndexMap map, double precursorMz, float fragmentBinSize, int length) {
		this.precursorMz=precursorMz;
		this.fragmentBinSize=fragmentBinSize;
		this.length=length;
		
		final ArrayList<SortablePeak> peaks=new ArrayList<SparseXCorrSpectrum.SortablePeak>();
		map.forEachEntry(new TIntObjectProcedure() {
			@Override
			public boolean execute(int a, Object b) {
				Peak peak=(Peak)b;
				peaks.add(new SortablePeak(a, peak.mass, peak.intensity));
				return true;
			}
		});
		Collections.sort(peaks);
		
		indices=new int[peaks.size()];
		masses=new double[peaks.size()];
		intensities=new float[peaks.size()];
		for (int i=0; i<indices.length; i++) {
			SortablePeak peak=peaks.get(i);
			indices[i]=peak.index;
			masses[i]=peak.mass;
			intensities[i]=peak.intensity;
		}
	}
	
	public double getPrecursorMz() {
		return precursorMz;
	}
	
	public float getFragmentBinSize() {
		return fragmentBinSize;
	}
	
	public int[] getIndices() {
		return indices;
	}
	
	public float[] getIntensityArray() {
		return intensities;
	}
	
	public double[] getMassArray() {
		return masses;
	}
	
	public int length() {
		return length;
	}
	
	public float[] toArray() {
		float[] array=new float[length];
		for (int i=0; i<indices.length; i++) {
			array[indices[i]]=intensities[i];
		}
		return array;
	}
	
	class SortablePeak implements Comparable<SortablePeak> {
		private final int index;
		private final double mass;
		private final float intensity;
		public SortablePeak(int index, double mass, float intensity) {
			this.index=index;
			this.mass=mass;
			this.intensity=intensity;
		}
		
		@Override
		public int compareTo(SortablePeak o) {
			if (o==null) return 1;
			return Integer.compare(index, o.index);
		}
	}

	public float dotProduct(SparseXCorrSpectrum spectrum) {
		int i=0;
		int j=0;
		float dotProduct=0.0f;
		while (i<indices.length&&j<spectrum.indices.length) {
			if (indices[i]==spectrum.indices[j]) {
				dotProduct+=intensities[i]*spectrum.intensities[j];
				i++;
				j++;
			} else if (indices[i]>spectrum.indices[j]) {
				j++;
			} else {
				i++;
			}
		}
		return dotProduct;
	}

	public float dotProduct(SparseXCorrSpectrum spectrum, int offset) {
		int i=0;
		int j=0;
		float dotProduct=0.0f;
		while (i<indices.length&&j<spectrum.indices.length) {
			int spectrumIndex=spectrum.indices[j]+offset;
			if (indices[i]==spectrumIndex) {
				dotProduct+=intensities[i]*spectrum.intensities[j];
				i++;
				j++;
			} else if (indices[i]>spectrumIndex) {
				j++;
			} else {
				i++;
			}
		}
		return dotProduct;
	}
}
