package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakChromatogram;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

//@Immutable
public class LibraryEntry implements Spectrum, PeptidePrecursor, XYTraceInterface {
	public static final String SHUFFLE_STRING="SHUFFLE_";
	public static final String DECOY_STRING="DECOY_";

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
	private final float[] correlationArray;
	private final HashSet<String> accessions;

	public LibraryEntry(String source, HashSet<String> accessions, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray) {
		this(source, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, getUnitArray(massArray.length));
	}
	
	private static float[] getUnitArray(int length) {
		float[] unit=new float[length];
		Arrays.fill(unit, 1.0f);
		return unit;
	}

	public LibraryEntry(String source, HashSet<String> accessions, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray, float[] correlationArray) {
		this(source, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, correlationArray);
	}

	public LibraryEntry(String source, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray) {
		this(source, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, getUnitArray(massArray.length));
		
	}
	public LibraryEntry(String source, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray, float[] correlationArray) {
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
			this.retentionTime=0.0f; //(float)SSRCalc.getHydrophobicity(peptideModSeq);
		}
		this.score=score;
		this.massArray=massArray;
		this.intensityArray=intensityArray;
		this.correlationArray=correlationArray;
	}
	
	@Override
	public Optional<Color> getColor() {
		return Optional.ofNullable((Color)null);
	}
	@Override
	public String getName() {
		return peptideModSeq;
	}
	@Override
	public Optional<Float> getThickness() {
		return Optional.ofNullable((Float)null);
	}
	@Override
	public GraphType getType() {
		return GraphType.spectrum;
	}
	@Override
	public Pair<double[], double[]> toArrays() {
		return new Pair<double[], double[]>(massArray, General.toDoubleArray(intensityArray));
	}
	
	public PercolatorPeptide getPSMData() {
		return new PercolatorPeptide(PercolatorPeptide.getPSMID(this, getRetentionTime(), new File(source)), PSMData.accessionsToString(accessions), getScore(), getScore());
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
	public LibraryEntry toUnitSpectrum(float rt) {
		return toUnitSpectrum(-1, rt);
	}
	public LibraryEntry toUnitSpectrum(int numPeaks) {
		return toUnitSpectrum(numPeaks, retentionTime);
	}
	public LibraryEntry toUnitSpectrum(int numPeaks, float rt) {
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
		return new LibraryEntry(source, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, rt, score, massArray, unit, correlationArray);
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
	public int compareTo(PeptidePrecursor o) {
		if (o==null) return 1;
		int c=peptideModSeq.compareTo(o.getPeptideModSeq());
		if (c!=0) return c;
		c=Byte.compare(precursorCharge, o.getPrecursorCharge());
		return c;
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

	public float[] getCorrelationArray() {
		return correlationArray;
	}
	
	public ArrayList<Peak> getPeaks() {
		ArrayList<Peak> peaks=new ArrayList<Peak>();
		for (int i = 0; i < massArray.length; i++) {
			peaks.add(new Peak(massArray[i], intensityArray[i]));
		}
		return peaks;
	}

	public LibraryEntry getDecoy(SearchParameters parameters) {
		return getDecoy(parameters, 0, false, true);
	} 
	public LibraryEntry getShuffle(SearchParameters parameters, int shuffleSeed, boolean markAsDecoy) {
		return getDecoy(parameters, shuffleSeed, true, markAsDecoy);
	} 
	private LibraryEntry getDecoy(SearchParameters parameters, int shuffleSeed, boolean shuffle, boolean markAsDecoy) {
		String reverseSequence;
		if (shuffle) {
			reverseSequence=PeptideUtils.shuffle(peptideModSeq, shuffleSeed, parameters);
		} else {
			reverseSequence=PeptideUtils.reverse(peptideModSeq, parameters);
		}
		
		FragmentationModel forwardModel=new FragmentationModel(peptideModSeq, parameters.getAAConstants());
		FragmentationModel reverseModel=new FragmentationModel(reverseSequence, parameters.getAAConstants());
		
		ArrayList<FragmentIon> forwardIons=new ArrayList<FragmentIon>();
		ArrayList<FragmentIon> reverseIons=new ArrayList<FragmentIon>();
		switch (parameters.getFragType()) {
		case YONLY:
			Collections.addAll(forwardIons, forwardModel.getYIons());
			Collections.addAll(reverseIons, reverseModel.getYIons());
			break;

		case CID:
			Collections.addAll(forwardIons, forwardModel.getBIons());
			Collections.addAll(reverseIons, reverseModel.getBIons());
			Collections.addAll(forwardIons, forwardModel.getYIons());
			Collections.addAll(reverseIons, reverseModel.getYIons());
			break;

		case ETD:
			Collections.addAll(forwardIons, forwardModel.getCIons());
			Collections.addAll(reverseIons, reverseModel.getCIons());
			Collections.addAll(forwardIons, forwardModel.getZIons());
			Collections.addAll(reverseIons, reverseModel.getZIons());
			Collections.addAll(forwardIons, forwardModel.getZp1Ions());
			Collections.addAll(reverseIons, reverseModel.getZp1Ions());
			break;
			
		}
		
		if (precursorCharge>2) {
			Collections.addAll(forwardIons, FragmentationModel.getPlus2s(forwardIons.toArray(new FragmentIon[forwardIons.size()])));
			Collections.addAll(reverseIons, FragmentationModel.getPlus2s(reverseIons.toArray(new FragmentIon[reverseIons.size()])));
		}
		
		assert(forwardIons.size()==reverseIons.size());
		ArrayList<XYPoint> points=new ArrayList<XYPoint>();
		int size=Math.min(forwardIons.size(), reverseIons.size()); //FIXME !!! HOW TO DEAL WITH LINKING UP NEUTRAL LOSSES IN DECOY PEPTIDES????
		for (int i=0; i<size; i++) {
			try {
				points.add(new XYPoint(forwardIons.get(i).mass, reverseIons.get(i).mass));
			} catch (Exception e) {
				System.out.println("WTF: "+peptideModSeq+"+"+precursorCharge+" vs "+reverseSequence+"+"+precursorCharge+", "+forwardIons.size()+" = "+forwardModel.getBIons().length+" + "+forwardModel.getYIons().length+"\t"+reverseIons.size()+" = "+reverseModel.getBIons().length+" + "+reverseModel.getYIons().length);
			
				throw new RuntimeException(e);
			}
		}
		Collections.sort(points);
		Pair<double[], double[]> matchedMasses=XYTrace.toArrays(points);
		double[] modelMasses=matchedMasses.x;
		double[] shiftedMasses=matchedMasses.y;
		
		float[] correlationArray=this.getCorrelationArray();

		MassTolerance tolerance=parameters.getFragmentTolerance();
		ArrayList<PeakChromatogram> reversedPeaks=new ArrayList<PeakChromatogram>();
		for (int i=0; i<massArray.length; i++) {
			double mass=massArray[i];
			float intensity=intensityArray[i];
			float correlation=correlationArray[i];
			
			Optional<Integer> matchIndex=tolerance.getIndex(modelMasses, mass);
			if (matchIndex.isPresent()) {
				double shiftedMass=shiftedMasses[matchIndex.get()];
				double delta=modelMasses[matchIndex.get()]-mass; // add back error if there is any
				
				// shift sequence specific ions
				reversedPeaks.add(new PeakChromatogram(shiftedMass-delta, intensity, correlation));
			} else {
				// add unknown peak with no modifications
				reversedPeaks.add(new PeakChromatogram(mass, intensity, correlation));
			}
		}
		Collections.sort(reversedPeaks);
		Triplet<double[], float[], float[]> arrays=PeakChromatogram.toChromatogramArrays(reversedPeaks);
		
		HashSet<String> revAcc=new HashSet<String>();
		for (String accession : accessions) {
			if (shuffle) {
				revAcc.add(DECOY_STRING+accession);
			} else {
				revAcc.add(SHUFFLE_STRING+accession);
			}
		}
		if (markAsDecoy) {
			return new ReverseLibraryEntry(source, revAcc, precursorMZ, precursorCharge, reverseSequence, copies, retentionTime, score, arrays.x, arrays.y, arrays.z);	
		} else {
			return new LibraryEntry(source, revAcc, precursorMZ, precursorCharge, reverseSequence, copies, retentionTime, score, arrays.x, arrays.y, arrays.z);	
		}
	}
}
