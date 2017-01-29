package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.Collections;

import gnu.trove.map.hash.TIntFloatHashMap;
import gnu.trove.procedure.TIntFloatProcedure;

public class SparseXCorrSpectrum {
	private final int[] indices;
	private final float[] intensities;
	private final int length;
	
	SparseXCorrSpectrum(TIntFloatHashMap map, int length) {
		this.length=length;
		
		final ArrayList<SortablePeak> peaks=new ArrayList<SparseXCorrSpectrum.SortablePeak>();
		map.forEachEntry(new TIntFloatProcedure() {	
			@Override
			public boolean execute(int a, float b) {
				peaks.add(new SortablePeak(a, b));
				return true;
			}
		});
		Collections.sort(peaks);
		
		indices=new int[peaks.size()];
		intensities=new float[peaks.size()];
		for (int i=0; i<indices.length; i++) {
			SortablePeak peak=peaks.get(i);
			indices[i]=peak.index;
			intensities[i]=peak.intensity;
		}
	}
	
	public int[] getIndices() {
		return indices;
	}
	
	public float[] getIntensities() {
		return intensities;
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
		private final float intensity;
		public SortablePeak(int index, float intensity) {
			this.index=index;
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
