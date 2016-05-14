package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.SSRCalc;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.list.array.TDoubleArrayList;

//@Immutable
public class LibraryEntry implements Comparable<LibraryEntry>, Spectrum {
	private static final float minimumIntensityThreshold=10.0f*Float.MIN_VALUE;
	
	private final String source;
	private final int spectrumIndex;
	private final double precursorMZ;
	private final byte precursorCharge;
	private final String peptideModSeq;
	private final int copies;
	private final float retentionTime;
	private final float score;
	private final double[] massArray;
	private final float[] intensityArray;
	private final HashSet<String> accessions;

	public LibraryEntry(String source, HashSet<String> accessions, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray) {
		this(source, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray);
	}

	public LibraryEntry(String source, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray) {
		this.source=source;
		this.accessions=new HashSet<String>(accessions);
		this.spectrumIndex=spectrumIndex;
		this.precursorMZ=precursorMZ;
		this.precursorCharge=precursorCharge;
		this.peptideModSeq=peptideModSeq;
		this.copies=copies;
		if (retentionTime>0.0f) { 
			this.retentionTime=retentionTime;
		} else {
			// ignores mods! This would be a problem if everything is modified (IMac prep)
			this.retentionTime=(float)SSRCalc.getHydrophobicity(peptideModSeq);
		}
		this.score=score;
		this.massArray=massArray;
		this.intensityArray=intensityArray;
	}
	
	public String getSource() {
		return source;
	}
	
	public HashSet<String> getAccessions() {
		return accessions;
	}

	public LibraryEntry toUnitSpectrum() {
		return toUnitSpectrum(-1);
	}
	public LibraryEntry toUnitSpectrum(int numPeaks) {
		float threshold;
		if (numPeaks<=0) {
			threshold=minimumIntensityThreshold;
		} else {
			float[] intensityArrayClone=intensityArray.clone();
			Arrays.sort(intensityArrayClone);
			int i=intensityArrayClone.length-numPeaks;
			if (i<0) i=0;
			threshold=intensityArrayClone[i];
		}
		
		float[] unit=new float[intensityArray.length];
		for (int i=0; i<unit.length; i++) {
			if (intensityArray[i]>=threshold) {
				unit[i]=1.0f;
			}
		}
		return new LibraryEntry(source, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, unit);
	}
	
	public float getTIC() {
		float tic=0.0f;
		for (int i=0; i<intensityArray.length; i++) {
			if (intensityArray[i]>minimumIntensityThreshold) {
				tic+=intensityArray[i];
			}
		}
		return tic;
	}
	
	public int getIonCount() {
		int count=0;
		for (int i=0; i<intensityArray.length; i++) {
			if (intensityArray[i]>minimumIntensityThreshold) {
				count++;
			}
		}
		return count;
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
	
	@Override
	public float getScanStartTime() {
		return getRetentionTime();
	}
	
	@Override
	public String getSpectrumName() {
		return peptideModSeq+"+"+precursorCharge;
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

	public ReverseLibraryEntry getDecoy(SearchParameters parameters, boolean shuffle) {
		String reverseSequence;
		if (shuffle) {
			reverseSequence=PeptideUtils.shuffle(peptideModSeq, parameters);
		} else {
			reverseSequence=PeptideUtils.reverse(peptideModSeq, parameters);
		}
		
		FragmentationModel forwardModel=new FragmentationModel(peptideModSeq, parameters.getAAConstants());
		FragmentationModel reverseModel=new FragmentationModel(reverseSequence, parameters.getAAConstants());
		
		TDoubleArrayList forwardIons=new TDoubleArrayList();
		TDoubleArrayList reverseIons=new TDoubleArrayList();
		switch (parameters.getFragType()) {
		case YONLY:
			forwardIons.add(forwardModel.getYIons());
			reverseIons.add(reverseModel.getYIons());
			break;

		case CID:
			forwardIons.add(forwardModel.getBIons());
			reverseIons.add(reverseModel.getBIons());
			forwardIons.add(forwardModel.getYIons());
			reverseIons.add(reverseModel.getYIons());
			break;

		case ETD:
			forwardIons.add(forwardModel.getCIons());
			reverseIons.add(reverseModel.getCIons());
			forwardIons.add(forwardModel.getZIons());
			reverseIons.add(reverseModel.getZIons());
			forwardIons.add(forwardModel.getZp1Ions());
			reverseIons.add(reverseModel.getZp1Ions());
			break;
			
		}
		
		if (precursorCharge>2) {
			forwardIons.add(FragmentationModel.getPlus2s(forwardIons.toArray()));
			reverseIons.add(FragmentationModel.getPlus2s(reverseIons.toArray()));
		}
		
		assert(forwardIons.size()==reverseIons.size());
		ArrayList<XYPoint> points=new ArrayList<XYPoint>();
		int size=Math.min(forwardIons.size(), reverseIons.size()); //FIXME !!! HOW TO DEAL WITH LINKING UP NEUTRAL LOSSES IN DECOY PEPTIES????
		for (int i=0; i<size; i++) {
			try {
			points.add(new XYPoint(forwardIons.get(i), reverseIons.get(i)));
			} catch (Exception e) {
				System.out.println("WTF: "+peptideModSeq+"+"+precursorCharge+" vs "+reverseSequence+"+"+precursorCharge+", "+forwardIons.size()+" = "+forwardModel.getBIons().length+" + "+forwardModel.getYIons().length+"\t"+reverseIons.size()+" = "+reverseModel.getBIons().length+" + "+reverseModel.getYIons().length);
			
				throw new RuntimeException(e);
			}
		}
		Collections.sort(points);
		Pair<double[], double[]> matchedMasses=XYTrace.toArrays(points);
		double[] modelMasses=matchedMasses.x;
		double[] shiftedMasses=matchedMasses.y;

		MassTolerance tolerance=parameters.getFragmentTolerance();
		ArrayList<Peak> reversedPeaks=new ArrayList<Peak>();
		for (int i=0; i<massArray.length; i++) {
			double mass=massArray[i];
			float intensity=intensityArray[i];
			
			Optional<Integer> matchIndex=tolerance.getIndex(modelMasses, mass);
			if (matchIndex.isPresent()) {
				double shiftedMass=shiftedMasses[matchIndex.get()];
				double delta=modelMasses[matchIndex.get()]-mass; // add back error if there is any
				
				// shift sequence specific ions
				reversedPeaks.add(new Peak(shiftedMass-delta, intensity));
			} else {
				// add unknown peak with no modifications
				reversedPeaks.add(new Peak(mass, intensity));
			}
		}
		Collections.sort(reversedPeaks);
		Pair<double[], float[]> arrays=Peak.toArrays(reversedPeaks);
		
		HashSet<String> revAcc=new HashSet<String>();
		for (String accession : accessions) {
			if (shuffle) {
				revAcc.add("DECOY_"+accession);
			} else {
				revAcc.add("SHUFFLE_"+accession);
			}
		}
		return new ReverseLibraryEntry(source, revAcc, precursorMZ, precursorCharge, reverseSequence, copies, retentionTime, score, arrays.x, arrays.y);	
	}
}
