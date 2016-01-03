package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Collections;

import com.google.common.base.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;

//@Immutable
public class LibraryEntry implements Comparable<LibraryEntry> {
	private final int spectrumIndex;
	private final double precursorMZ;
	private final byte precursorCharge;
	private final String peptideModSeq;
	private final int copies;
	private final float retentionTime;
	private final float score;
	private final double[] massArray;
	private final float[] intensityArray;


	public LibraryEntry(double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray) {
		this(1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray);
	}

	public LibraryEntry(int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray) {
		this.spectrumIndex=spectrumIndex;
		this.precursorMZ=precursorMZ;
		this.precursorCharge=precursorCharge;
		this.peptideModSeq=peptideModSeq;
		this.copies=copies;
		this.retentionTime=retentionTime;
		this.score=score;
		this.massArray=massArray;
		this.intensityArray=intensityArray;
	}
	
	@Override
	public int compareTo(LibraryEntry o) {
		if (o==null) return 1;
		int c=peptideModSeq.compareTo(o.peptideModSeq);
		if (c!=0) return c;
		c=Byte.compare(precursorCharge, o.precursorCharge);
		if (c!=0) return c;
		return Float.compare(retentionTime, o.retentionTime);
	}
	
	public int getSpectrumIndex() {
		return spectrumIndex;
	}
	
	public boolean isDecoy() {
		return false;
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
	
	public String getPeptideSeq() {
		StringBuilder sb=new StringBuilder();
		for (char c : peptideModSeq.toCharArray()) {
			if (Character.isLetter(c)) {
				sb.append(c);
			}
		}
		return sb.toString();
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
	
	public ReverseLibraryEntry getReverse(MassTolerance tolerance, AminoAcidConstants aaConstants) {
		FragmentationModel model=new FragmentationModel(peptideModSeq, aaConstants);
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
