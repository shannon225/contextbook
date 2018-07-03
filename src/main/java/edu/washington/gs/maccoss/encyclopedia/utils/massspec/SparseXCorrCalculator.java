package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.SparseIndexMap;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class SparseXCorrCalculator {
	public static final double biggestFragmentMass=2000.0;
	
	private final SearchParameters params;
	private final SparseXCorrSpectrum preprocessedSpectrum;
	
	/**
	 * 
	 * @param s assumes sqrted intensities!
	 * @param precursorMz
	 * @param charge
	 * @param params
	 */
	public SparseXCorrCalculator(Spectrum s, Range precursorMz, SearchParameters params) {
		this(normalize(s, precursorMz, false, params), params);
	}
	
	public SparseXCorrCalculator(String modifiedSequence, byte precursorCharge, SearchParameters params) {
		this(getTheoreticalSpectrum(modifiedSequence, precursorCharge, params), params);
	}
	
	public SparseXCorrCalculator(SparseXCorrSpectrum spectrum, SearchParameters params) {
		this.params=params;
		preprocessedSpectrum=preprocessSpectrum(spectrum);
	}
	
	/**
	 * 
	 * @param s assumes sqrted intensities!
	 * @return
	 */
	public float score(Spectrum s, Range precursorMz) {
		SparseXCorrSpectrum intensityBins=normalize(s, precursorMz);
		return score(intensityBins);
	}

	public SparseXCorrSpectrum normalize(Spectrum s, Range precursorMz) {
		return normalize(s, precursorMz, false, params);
	}
	
	public float score(String modifiedSequence, byte precursorCharge) {
		SparseXCorrSpectrum intensityBins=getTheoreticalSpectrum(modifiedSequence, precursorCharge, params);
		return score(intensityBins);
	}
	
	public float score(SparseXCorrSpectrum spectrum) {
		// divide by 1e4 (personal communication with J Eng)
		return preprocessedSpectrum.dotProduct(spectrum)/1.0e4f;
	}

	static float dotProduct(float[] preprocessedSpectrum, float[] spectrum) {
		float sum=0.0f;
		for (int i=0; i<spectrum.length; i++) {
			if (i>=preprocessedSpectrum.length) break;
			sum+=spectrum[i]*preprocessedSpectrum[i];
		}
		return sum;
	}

	static float dotProduct(float[] preprocessedSpectrum, float[] spectrum, int offset) {
		float sum=0.0f;
		for (int i=0; i<spectrum.length; i++) {
			int index=i+offset;
			if (index<0||index>=preprocessedSpectrum.length) continue;
			sum+=spectrum[i]*preprocessedSpectrum[index];
		}
		return sum;
	}

	static SparseXCorrSpectrum preprocessSpectrum(SparseXCorrSpectrum spectrum) {
		SparseIndexMap preprocessedSpectrum=new SparseIndexMap();
		
		int length=spectrum.length();
		int[] indicies=spectrum.getIndices();
		double[] masses=spectrum.getMassArray();
		float[] intensities=spectrum.getIntensityArray();
		float[] negativeIntensities=new float[intensities.length];
		for (int i=0; i<negativeIntensities.length; i++) {
			negativeIntensities[i]=-intensities[i];
		}
		
		for (int offset=ArrayXCorrCalculator.lowerOffset; offset<ArrayXCorrCalculator.upperOffset; offset++) {
			if (offset==0) continue;
			for (int i=0; i<indicies.length; i++) {
				int index=indicies[i]+offset;

				if (index>=0&&index<length) {
					preprocessedSpectrum.adjustOrPutValue(index, masses[i]+offset*spectrum.getFragmentBinSize(), negativeIntensities[i]);
				}
			}
		}

		final int denominator=ArrayXCorrCalculator.upperOffset-ArrayXCorrCalculator.lowerOffset;
		preprocessedSpectrum.multiplyAllValues(1.0f/denominator);
		
		for (int i=0; i<indicies.length; i++) {
			preprocessedSpectrum.adjustOrPutValue(indicies[i], masses[i], intensities[i]);
		}
		
		return new SparseXCorrSpectrum(preprocessedSpectrum, spectrum.getPrecursorMZ(), spectrum.getFragmentBinSize(), length);
	}
	
	/**
	 * see Eng et al, JASMS 1994
	 * @param s
	 * @param precursorMz
	 * @return
	 */
	public static SparseXCorrSpectrum normalize(Spectrum s, Range precursorMz, boolean addIntensityToNeighboringBins, SearchParameters params) {
		
		double[] masses=s.getMassArray();
		float[] intensities=s.getIntensityArray();
		ArrayList<Peak> allPeaks=new ArrayList<Peak>();
		if (masses.length==0)
			return getIntensityArray(params, allPeaks, s.getPrecursorMZ(), addIntensityToNeighboringBins);
		if (masses.length==1) {
			allPeaks.add(new Peak(masses[0], ArrayXCorrCalculator.primaryIonIntensity));
			return getIntensityArray(params, allPeaks, s.getPrecursorMZ(), addIntensityToNeighboringBins);
		}

		double minimumPrecursorRemoved=precursorMz.getStart();
		double maximumPrecursorRemoved=precursorMz.getStop();

		double firstMass=masses[0];
		double lastMass=masses[masses.length-1];

		double increment=(lastMass-firstMass)/ArrayXCorrCalculator.groups;
		double[] binMaxMass=new double[ArrayXCorrCalculator.groups]; 
		for (int i=0; i<ArrayXCorrCalculator.groups-1; i++) {
			binMaxMass[i]=increment*(i+1);
		}
		binMaxMass[ArrayXCorrCalculator.groups-1]=Double.MAX_VALUE;
		
		float[] binMaxIntensity=new float[ArrayXCorrCalculator.groups];
		int currentIndex=0;
		for (int i=0; i<intensities.length; i++) {
			if (intensities[i]<=0.0f||masses[i]>minimumPrecursorRemoved&&masses[i]<maximumPrecursorRemoved) {
				continue;
			}
			
			while (masses[i]>binMaxMass[currentIndex]) {
				currentIndex++;
			}
			
			if (intensities[i]>binMaxIntensity[currentIndex]) {
				binMaxIntensity[currentIndex]=intensities[i];
			}
		}
		
		binMaxIntensity=General.divide(binMaxIntensity, ArrayXCorrCalculator.primaryIonIntensity);
		
		currentIndex=0;
		for (int i=0; i<intensities.length; i++) {
			if (intensities[i]<=0.0f||(masses[i]>minimumPrecursorRemoved&&masses[i]<maximumPrecursorRemoved)) {
				continue;
			}
			
			while (masses[i]>binMaxMass[currentIndex]) {
				currentIndex++;
			}
			allPeaks.add(new Peak(masses[i], intensities[i]/binMaxIntensity[currentIndex]));
		}
		
		return getIntensityArray(params, allPeaks, s.getPrecursorMZ(), addIntensityToNeighboringBins);
	}
	
	public static SparseXCorrSpectrum getTheoreticalSpectrum(String modifiedSequence, byte precursorCharge, SearchParameters params) {
		return getTheoreticalSpectrumPair(modifiedSequence, precursorCharge, params).y;
	}

	public static Pair<FragmentationModel, SparseXCorrSpectrum> getTheoreticalSpectrumPair(String modifiedSequence, byte precursorCharge, SearchParameters params) {
		
		FragmentationType type=params.getFragType();
		AminoAcidConstants aaConstants=params.getAAConstants();
		FragmentationModel model=PeptideUtils.getPeptideModel(modifiedSequence, aaConstants);
		
		ArrayList<Peak> allPeaks=new ArrayList<Peak>();
		switch (type) {
			case HCD:
				FragmentIon[] yIons=model.getYIons();
				allPeaks.addAll(getPeaks(yIons, 0.0, ArrayXCorrCalculator.primaryIonIntensity));
				if (params.isUseNLsForXCorr()) {
					allPeaks.addAll(getPeaks(yIons, -MassConstants.nh3, ArrayXCorrCalculator.neutralLossIntensity));
					allPeaks.addAll(getPeaks(yIons, -MassConstants.oh2, ArrayXCorrCalculator.neutralLossIntensity));
				}
				break;
				
			case CID:
				FragmentIon[] yIonsCID=model.getYIons();
				allPeaks.addAll(getPeaks(yIonsCID, 0.0, ArrayXCorrCalculator.primaryIonIntensity));
				if (params.isUseNLsForXCorr()) {
					allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.nh3, ArrayXCorrCalculator.neutralLossIntensity));
					allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.oh2, ArrayXCorrCalculator.neutralLossIntensity));
				}
				
				FragmentIon[] bIonsCID=model.getBIons();
				allPeaks.addAll(getPeaks(bIonsCID, 0.0, ArrayXCorrCalculator.primaryIonIntensity));
				if (params.isUseNLsForXCorr()) {
					allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.nh3, ArrayXCorrCalculator.neutralLossIntensity));
					allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.oh2, ArrayXCorrCalculator.neutralLossIntensity));
					allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.co, ArrayXCorrCalculator.neutralLossIntensity));
				}
				break;
				
			case ETD:
				FragmentIon[] cIonsCID=model.getCIons();
				allPeaks.addAll(getPeaks(cIonsCID, 0.0, ArrayXCorrCalculator.primaryIonIntensity));
				if (params.isUseNLsForXCorr()) {
					allPeaks.addAll(getPeaks(cIonsCID, -MassConstants.nh3, ArrayXCorrCalculator.neutralLossIntensity));
					allPeaks.addAll(getPeaks(cIonsCID, -MassConstants.oh2, ArrayXCorrCalculator.neutralLossIntensity));
				}

				FragmentIon[] zIonsCID=model.getCIons();
				allPeaks.addAll(getPeaks(zIonsCID, 0.0, ArrayXCorrCalculator.primaryIonIntensity));
				allPeaks.addAll(getPeaks(zIonsCID, MassConstants.neutronMass, ArrayXCorrCalculator.primaryIonIntensity)); // z+1
				if (params.isUseNLsForXCorr()) {
					allPeaks.addAll(getPeaks(zIonsCID, -MassConstants.nh3, ArrayXCorrCalculator.neutralLossIntensity));
					allPeaks.addAll(getPeaks(zIonsCID, -MassConstants.oh2, ArrayXCorrCalculator.neutralLossIntensity));
				}
				break;
				
			default:
				throw new EncyclopediaException("Unknown fragmentation type ["+type+"]");
		}
		
		return new Pair<>(model, getIntensityArray(params, allPeaks, model.getChargedMass(precursorCharge), true));
	}
	
	private static ArrayList<Peak> getPeaks(FragmentIon[] ions, double delta, float intensity) {
		ArrayList<Peak> peaks=new ArrayList<Peak>();
		for (int i=0; i<ions.length; i++) {
			peaks.add(new Peak(ions[i].mass+delta, intensity));
		}
		return peaks;
	}

	private static SparseXCorrSpectrum getIntensityArray(SearchParameters params, ArrayList<Peak> allPeaks, double precursorMz, boolean addIntensityToNeighboringBins) {
		Collections.sort(allPeaks);
		
		// set tolerance to 2x the fragment tolerance of the highest fragment
		float fragmentBinSize=2.0f*(float)params.getFragmentTolerance().getTolerance(biggestFragmentMass);
		double offset;
		
		if (fragmentBinSize>0.5f) {
			fragmentBinSize=ArrayXCorrCalculator.lowResFragmentBinSize; // if tolerance is >0.25 Da, then jump to 1 Da to make use of the average amino acid mass defect
			offset=ArrayXCorrCalculator.lowResFragmentBinOffset;
		} else if (fragmentBinSize<0.01f) {
			fragmentBinSize=0.01f;
			offset=0.0;
		} else {
			offset=0.0;
		}

		float inverseBinWidth=1.0f/fragmentBinSize;
		int arraySize=(int)((biggestFragmentMass+fragmentBinSize+2.0)*inverseBinWidth);
		
		SparseIndexMap binnedIntensityArray=new SparseIndexMap(allPeaks.size());
		int arraySizeMinusOne=arraySize-1;
		for (Peak peak : allPeaks) {
			int massIndex=(int)((peak.mass-offset)*inverseBinWidth);
			
			if (massIndex<0) massIndex=0;
			if (massIndex>=arraySize) massIndex=arraySize-1;
			
			binnedIntensityArray.putIfGreater(massIndex, peak.mass, peak.intensity);
			
			// don't do this for low res fragment ions bin boundaries aren't an issue with the 0.4 offset
			if (fragmentBinSize<=0.5f&&addIntensityToNeighboringBins) {
				// neighboring intensities are 25 for b/y or 10 (the same) for neutral losses
				float neighboringIntensity=peak.intensity>ArrayXCorrCalculator.neutralLossIntensity?peak.intensity/2.0f:peak.intensity;
				if (massIndex>0) {
					binnedIntensityArray.putIfGreater(massIndex-1, peak.mass-fragmentBinSize, neighboringIntensity);
				}
				if (massIndex<arraySizeMinusOne) {
					binnedIntensityArray.putIfGreater(massIndex+1, peak.mass+fragmentBinSize, neighboringIntensity);
				}
			}
		}
		return new SparseXCorrSpectrum(binnedIntensityArray, precursorMz, fragmentBinSize, arraySize);
	}
}
