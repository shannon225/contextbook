package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IndexedIonType;
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
	private final String massCorrectedPeptideModSeq;
	private final int copies;
	private final float retentionTime;
	private final float score;
	private final double[] massArray;
	private final float[] intensityArray;
	private final float[] correlationArray;
	private final HashSet<String> accessions;

	public LibraryEntry(String source, HashSet<String> accessions, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray, AminoAcidConstants aaConstants) {
		this(source, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, getUnitArray(massArray.length), aaConstants);
	}
	
	private static float[] getUnitArray(int length) {
		float[] unit=new float[length];
		Arrays.fill(unit, 1.0f);
		return unit;
	}

	public LibraryEntry(String source, HashSet<String> accessions, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray, float[] correlationArray, AminoAcidConstants aaConstants) {
		this(source, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, correlationArray, aaConstants);
	}

	public LibraryEntry(String source, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray, AminoAcidConstants aaConstants) {
		this(source, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, getUnitArray(massArray.length), aaConstants);
	}

	public LibraryEntry(String source, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray, float[] correlationArray, AminoAcidConstants aaConstants) {
		this(
				source,
				accessions,
				spectrumIndex,
				precursorMZ,
				precursorCharge,
				peptideModSeq,
				PeptideUtils.getCorrectedMasses(peptideModSeq, aaConstants),
				copies,
				retentionTime,
				score,
				massArray,
				intensityArray,
				correlationArray
		);
	}

	public LibraryEntry(String source, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, String massCorrectedPeptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray, float[] correlationArray) {
		this.source=source;
		this.accessions=new HashSet<String>(accessions);
		this.spectrumIndex=spectrumIndex;
		this.precursorMZ=precursorMZ;
		this.precursorCharge=precursorCharge;
		this.peptideModSeq=peptideModSeq;
		this.massCorrectedPeptideModSeq=massCorrectedPeptideModSeq;
		this.copies=copies;
		if (retentionTime!=0.0f) { 
			this.retentionTime=retentionTime;
		} else {
			// ignores mods! This would be a problem if everything is modified (IMac prep)
			this.retentionTime=0.0f; //(float)SSRCalc.getHydrophobicity(peptideModSeq);
		}
		this.score=score;
		
		ArrayList<PeakChromatogram> peaks=new ArrayList<>();
		int numPeaks=Math.min(massArray.length, correlationArray.length);
		for (int i=0; i<numPeaks; i++) {
			peaks.add(new PeakChromatogram(massArray[i], intensityArray[i], correlationArray[i]));
		}
		Collections.sort(peaks);
		Triplet<double[], float[], float[]> arrays=PeakChromatogram.toChromatogramArrays(peaks);
		
		this.massArray=arrays.x;
		this.intensityArray=arrays.y;
		this.correlationArray=arrays.z;
	}
	
	@Override
	public String getPeptideModSeq() {
		return massCorrectedPeptideModSeq;
	}
	
	/**
	 * only use for testing
	 * @param rtInSec
	 * @return
	 */
	public LibraryEntry updateRetentionTime(float rtInSec) {
		return new LibraryEntry(source, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, massCorrectedPeptideModSeq, copies, rtInSec, score, massArray, intensityArray, correlationArray);
	}
	
	/**
	 * only use for testing
	 * @param newMassArray
	 * @param newIntensityArray
	 * @return
	 */
	public LibraryEntry updateMS2(double[] newMassArray, float[] newIntensityArray) {
		return new LibraryEntry(source, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, massCorrectedPeptideModSeq, copies, retentionTime, score, newMassArray, newIntensityArray, getUnitArray(newMassArray.length));
	}
	
	/**
	 * only use for testing
	 */
	public String toObjectCreatorString() {
		StringBuilder sb=new StringBuilder();
		sb.append("double[] massArray=new double[] {");
		for (int i=0; i<massArray.length; i++) {
			if (i>0) sb.append(", ");
			sb.append(massArray[i]);
		}
		sb.append("};");
		sb.append("float[] intensityArray=new float[] {");
		for (int i=0; i<intensityArray.length; i++) {
			if (i>0) sb.append(", ");
			sb.append(intensityArray[i]+"f");
		}
		sb.append("};");
		sb.append("float[] correlationArray=new float[] {");
		for (int i=0; i<correlationArray.length; i++) {
			if (i>0) sb.append(", ");
			sb.append(correlationArray[i]+"f");
		}
		sb.append("};\n");
		
		sb.append("LibraryEntry "+getPeptideSeq()+"=new LibraryEntry(\""+source+"\", new HashSet<String>(), "+spectrumIndex+", "+precursorMZ+", (byte)"+precursorCharge+", \""+peptideModSeq+"\", "+copies+", "+retentionTime+"f, "+score+"f, massArray, intensityArray, correlationArray, aaConstants);");
		return sb.toString();
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
	
	@Override
	public int size() {
		return massArray.length;
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
		return new LibraryEntry(source, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, massCorrectedPeptideModSeq, copies, rt, score, massArray, unit, correlationArray);
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
		int c=getPeptideModSeq().compareTo(o.getPeptideModSeq());
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

	public String getLegacyPeptideModSeq() {
		return peptideModSeq;
	}
	
	volatile String accuratePeptideModSeq=null;
	/**
	 * safe in that it will always return not null but it may recalculate accuratePeptideModSeq if there is thread competition
	 * @param aaConstants
	 * @return
	 */
	public String getAccuratePeptideModSeq(AminoAcidConstants aaConstants) {
		String sequence=accuratePeptideModSeq;
		if (sequence==null) {
			sequence=PeptideUtils.getSequence(PeptideUtils.getPeptideModel(peptideModSeq, aaConstants).getAas());
			accuratePeptideModSeq=sequence;
			return sequence;
		} 
		return sequence;
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
			reverseSequence=PeptideUtils.reverse(peptideModSeq, parameters.getAAConstants());
		}
		HashSet<String> revAcc=new HashSet<String>();
		for (String accession : accessions) {
			if (shuffle) {
				revAcc.add(SHUFFLE_STRING+accession);
			} else {
				revAcc.add(DECOY_STRING+accession);
			}
		}
		
		return getEntryFromNewSequence(reverseSequence, revAcc, markAsDecoy, parameters).y;
	}

	public Pair<FragmentationModel, LibraryEntry> getEntryFromNewSequence(String newSequence, HashSet<String> accessions, boolean markAsDecoy, SearchParameters parameters) {
		FragmentationModel forwardModel=PeptideUtils.getPeptideModel(peptideModSeq, parameters.getAAConstants());
		FragmentationModel reverseModel=PeptideUtils.getPeptideModel(newSequence, parameters.getAAConstants());

		double oldPrecursorMz=forwardModel.getChargedMass(precursorCharge);
		double precursorMzError=precursorMZ-oldPrecursorMz;
		double newPrecursorMz=reverseModel.getChargedMass(precursorCharge)-precursorMzError;
		
		ArrayList<FragmentIon> forwardIons=new ArrayList<FragmentIon>();
		ArrayList<FragmentIon> reverseIons=new ArrayList<FragmentIon>();
		switch (parameters.getFragType()) {
		case HCD:
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
		
		HashMap<IndexedIonType, FragmentIon> forwardMap=new HashMap<>();
		for (FragmentIon ion : forwardIons) {
			forwardMap.put(new IndexedIonType(ion), ion);
		}

		// make sure ion indices line up
		ArrayList<XYPoint> points=new ArrayList<XYPoint>();
		for (FragmentIon reverse : reverseIons) {
			FragmentIon forward=forwardMap.get(new IndexedIonType(reverse));
			if (forward!=null) {
				points.add(new XYPoint(forward.getMass(), reverse.getMass()));
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
		
		if (markAsDecoy) {
			return new Pair<FragmentationModel, LibraryEntry>(reverseModel, new ReverseLibraryEntry(source, accessions, newPrecursorMz, precursorCharge, newSequence, copies, retentionTime, score, arrays.x, arrays.y, arrays.z, parameters.getAAConstants()));
		} else {
			return new Pair<FragmentationModel, LibraryEntry>(reverseModel, new LibraryEntry(source, accessions, newPrecursorMz, precursorCharge, newSequence, copies, retentionTime, score, arrays.x, arrays.y, arrays.z, parameters.getAAConstants()));
		}
	}
}
