package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Collections;

import com.google.common.base.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import jdk.nashorn.internal.ir.annotations.Immutable;

@Immutable
public class LibraryEntry {

	private final double precursorMZ;
	private final byte precursorCharge;
	private final String peptideModSeq;
	private final int copies;
	private final float retentionTime;
	private final float score;
	private final double[] massArray;
	private final float[] intensityArray;

	public LibraryEntry(double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray) {
		this.precursorMZ=precursorMZ;
		this.precursorCharge=precursorCharge;
		this.peptideModSeq=peptideModSeq;
		this.copies=copies;
		this.retentionTime=retentionTime;
		this.score=score;
		this.massArray=massArray;
		this.intensityArray=intensityArray;
	}

	public double getPrecursorMZ() {
		return precursorMZ;
	}

	public byte getPrecursorCharge() {
		return precursorCharge;
	}

	public String getPeptideModSeq() {
		return peptideModSeq;
	}

	public int getCopies() {
		return copies;
	}

	public float getRetentionTime() {
		return retentionTime;
	}

	public float getScore() {
		return score;
	}

	public double[] getMassArray() {
		return massArray;
	}

	public float[] getIntensityArray() {
		return intensityArray;
	}
	
	public ArrayList<Peak> getPeaks() {
		ArrayList<Peak> peaks=new ArrayList<Peak>();
		for (int i = 0; i < massArray.length; i++) {
			peaks.add(new Peak(massArray[i], intensityArray[i]));
		}
		return peaks;
	}
	
	public ReverseLibraryEntry getReverse(MassTolerance tolerance) {
		FragmentationModel model=new FragmentationModel(peptideModSeq);
		double[] bs=model.getBIons();
		double[] ys=model.getYIons();
		
		ArrayList<Peak> reversedPeaks=new ArrayList<Peak>();
		for (int i = 0; i < massArray.length; i++) {
			// check if b ion
			Optional<Double> match = tolerance.getMatch(bs, massArray[i]);
			if (match.isPresent()) {
				// b+18.01057=reversed y
				reversedPeaks.add(new Peak(match.get()+18.01057, intensityArray[i]));
				continue;
			}
			
			// check if y ion
			match = tolerance.getMatch(ys, massArray[i]);
			if (match.isPresent()) {
				// y-18.01057=reversed b
				reversedPeaks.add(new Peak(match.get()-18.01057, intensityArray[i]));
				continue;
			}
			
			// if not b or y, then keep as is
			reversedPeaks.add(new Peak(massArray[i], intensityArray[i]));
		}
		
		Collections.sort(reversedPeaks);
		Pair<double[], float[]> arrays=Peak.toArrays(reversedPeaks);
		
		StringBuilder sb=new StringBuilder();
		String[] aas=model.getAas();
		for (int i = aas.length-1; i >=0; i--) {
			sb.append(aas[i]);
		}

		return new ReverseLibraryEntry(precursorMZ, precursorCharge, sb.toString(), copies, retentionTime, score, arrays.x, arrays.y);	
	}
}
