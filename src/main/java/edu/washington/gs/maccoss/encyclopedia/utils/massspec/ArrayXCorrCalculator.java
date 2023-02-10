package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class ArrayXCorrCalculator {
	// values from personal communication with J. Egertson
	static final float lowResFragmentBinSize=1.00045475f;
	static final float lowResFragmentBinOffset=0.4f;
	static final float highResFragmentBinSize=0.02f;
	static final float highResFragmentBinOffset=0.0f;
	
	// set 50 to be the maximum value, see pp 982 bottom right
	static float primaryIonIntensity=50.0f;
	static final float neutralLossIntensity=primaryIonIntensity/5;

	// offset defined figure legend on pp 980
	static final int upperOffset=75;
	static final int lowerOffset=-upperOffset;

	// divide spectrum into 10 equal regions, see pp 982 bottom right
	static int groups=10; 
	
	// remove 10-u window around precursor, see pp 979 mid left
	static double precursorRemovalMargin=5.0;
	
	private final SearchParameters params;
	private final byte charge;
	private final double precursorMz;
	private final float[] preprocessedSpectrum;
	
	/**
	 * 
	 * @param s assumes sqrted intensities!
	 * @param precursorMz
	 * @param charge
	 * @param params
	 */
	public ArrayXCorrCalculator(Spectrum s, double precursorMz, byte charge, SearchParameters params) {
		this.precursorMz=precursorMz;
		this.charge=charge;
		this.params=params;
		preprocessedSpectrum=preprocessSpectrum(normalize(s, precursorMz, charge, false, params));
	}
	
	public ArrayXCorrCalculator(String modifiedSequence, double precursorMz, byte charge, SearchParameters params) {
		this.precursorMz=precursorMz;
		this.charge=charge;
		this.params=params;
		preprocessedSpectrum=preprocessSpectrum(getTheoreticalSpectrum(modifiedSequence, precursorMz, charge, params));
	}
	
	/**
	 * 
	 * @param s assumes sqrted intensities!
	 * @return
	 */
	public float score(Spectrum s) {
		float[] intensityBins=normalize(s);
		return score(intensityBins);
	}

	public float[] normalize(Spectrum s) {
		return normalize(s, precursorMz, charge, false, params);
	}
	
	public float score(String modifiedSequence) {
		float[] intensityBins=getTheoreticalSpectrum(modifiedSequence, precursorMz, charge, params);
		return score(intensityBins);
	}
	
	/**
	 * divide by 1e4 (personal communication with J Eng)
	 * @param spectrum
	 * @return
	 */
	float score(float[] spectrum) {
		return dotProduct(preprocessedSpectrum, spectrum)/1.0e4f;
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

	static float[] preprocessSpectrum(float[] spectrum) {
		float[] preprocessedSpectrum=new float[spectrum.length];
		
		for (int offset=lowerOffset; offset<upperOffset; offset++) {
			if (offset==0) continue;
			for (int i=0; i<spectrum.length; i++) {
				int index=i+offset;
				if (index>=0&&index<preprocessedSpectrum.length) {
					preprocessedSpectrum[i]-=spectrum[index];
				}
			}
		}

		int denominator=upperOffset-lowerOffset;
		for (int i=0; i<preprocessedSpectrum.length; i++) {
			preprocessedSpectrum[i]=preprocessedSpectrum[i]/denominator;
		}
		
		for (int i=0; i<spectrum.length; i++) {
			preprocessedSpectrum[i]+=spectrum[i];
		}
		return preprocessedSpectrum;
	}
	
	/**
	 * see Eng et al, JASMS 1994
	 * @param s
	 * @param precursorMz
	 * @return
	 */
	static float[] normalize(Spectrum s, double precursorMz, byte charge, boolean addIntensityToNeighboringBins, SearchParameters params) {
		double massPlusOne=precursorMz*charge-(charge-1)*MassConstants.protonMass;
		
		double[] masses=s.getMassArray();
		float[] intensities=s.getIntensityArray();
		ArrayList<Peak> allPeaks=new ArrayList<Peak>();
		if (masses.length==0)
			return getIntensityArray(params, allPeaks, massPlusOne, addIntensityToNeighboringBins);
		if (masses.length==1) {
			allPeaks.add(new Peak(masses[0], primaryIonIntensity));
			return getIntensityArray(params, allPeaks, massPlusOne, addIntensityToNeighboringBins);
		}

		double minimumPrecursorRemoved=precursorMz-precursorRemovalMargin;
		double maximumPrecursorRemoved=precursorMz+precursorRemovalMargin;

		double firstMass=masses[0];
		double lastMass=masses[masses.length-1];

		double increment=(lastMass-firstMass)/groups;
		double[] binMaxMass=new double[groups]; 
		for (int i=0; i<groups-1; i++) {
			binMaxMass[i]=increment*(i+1);
		}
		binMaxMass[groups-1]=Double.MAX_VALUE;
		
		float[] binMaxIntensity=new float[groups];
		int currentIndex=0;
		for (int i=0; i<intensities.length; i++) {
			if (masses[i]>minimumPrecursorRemoved&&masses[i]<maximumPrecursorRemoved) {
				continue;
			}
			
			while (masses[i]>binMaxMass[currentIndex]) {
				currentIndex++;
			}
			
			if (intensities[i]>binMaxIntensity[currentIndex]) {
				binMaxIntensity[currentIndex]=intensities[i];
			}
		}
		
		binMaxIntensity=General.divide(binMaxIntensity, primaryIonIntensity);
		
		currentIndex=0;
		for (int i=0; i<intensities.length; i++) {
			if (masses[i]>minimumPrecursorRemoved&&masses[i]<maximumPrecursorRemoved) {
				continue;
			}
			
			while (masses[i]>binMaxMass[currentIndex]) {
				currentIndex++;
			}
			allPeaks.add(new Peak(masses[i], intensities[i]/binMaxIntensity[currentIndex]));
		}
		
		return getIntensityArray(params, allPeaks, massPlusOne, addIntensityToNeighboringBins);
	}
	
	static float[] getTheoreticalSpectrum(String modifiedSequence, double precursorMz, byte charge, SearchParameters params) {
		double massPlusOne=precursorMz*charge-(charge-1)*MassConstants.protonMass;
		
		FragmentationType type=params.getFragType();
		AminoAcidConstants aaConstants=params.getAAConstants();
		FragmentationModel model=PeptideUtils.getPeptideModel(modifiedSequence, aaConstants);
		
		ArrayList<Peak> allPeaks=new ArrayList<Peak>();
		switch (type) {
			case HCD:
			case SILAC:
				Ion[] yIons=model.getYIons();
				allPeaks.addAll(getPeaks(yIons, 0.0, primaryIonIntensity));
				if (params.isUseNLsForXCorr()) {
					allPeaks.addAll(getPeaks(yIons, -MassConstants.nh3, neutralLossIntensity));
					allPeaks.addAll(getPeaks(yIons, -MassConstants.oh2, neutralLossIntensity));
				}
				break;

			case CID:
				Ion[] yIonsCID=model.getYIons();
				allPeaks.addAll(getPeaks(yIonsCID, 0.0, primaryIonIntensity));
				if (params.isUseNLsForXCorr()) {
					allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.nh3, neutralLossIntensity));
					allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.oh2, neutralLossIntensity));
				}

				Ion[] bIonsCID=model.getBIons();
				allPeaks.addAll(getPeaks(bIonsCID, 0.0, primaryIonIntensity));
				if (params.isUseNLsForXCorr()) {
					allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.nh3, neutralLossIntensity));
					allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.oh2, neutralLossIntensity));
					allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.co, neutralLossIntensity));
				}
				break;

			case ETD:
				Ion[] cIonsCID=model.getCIons();
				allPeaks.addAll(getPeaks(cIonsCID, 0.0, primaryIonIntensity));
				if (params.isUseNLsForXCorr()) {
					allPeaks.addAll(getPeaks(cIonsCID, -MassConstants.nh3, neutralLossIntensity));
					allPeaks.addAll(getPeaks(cIonsCID, -MassConstants.oh2, neutralLossIntensity));
				}

				Ion[] zIonsCID=model.getCIons();
				allPeaks.addAll(getPeaks(zIonsCID, 0.0, primaryIonIntensity));
				allPeaks.addAll(getPeaks(zIonsCID, MassConstants.neutronMass, primaryIonIntensity)); // z+1
				if (params.isUseNLsForXCorr()) {
					allPeaks.addAll(getPeaks(zIonsCID, -MassConstants.nh3, neutralLossIntensity));
					allPeaks.addAll(getPeaks(zIonsCID, -MassConstants.oh2, neutralLossIntensity));
				}
				break;
				
			default:
				throw new EncyclopediaException("Unknown fragmentation type ["+type+"]");
		}
		
		return getIntensityArray(params, allPeaks, massPlusOne, true);
	}
	
	private static ArrayList<Peak> getPeaks(Ion[] ions, double delta, float intensity) {
		ArrayList<Peak> peaks=new ArrayList<Peak>();
		for (int i=0; i<ions.length; i++) {
			peaks.add(new Peak(ions[i].getMass()+delta, intensity));
		}
		return peaks;
	}

	private static float[] getIntensityArray(SearchParameters params, ArrayList<Peak> allPeaks, double massPlusOne, boolean addIntensityToNeighboringBins) {
		Collections.sort(allPeaks);
		
		// set tolerance to 2x the fragment tolerance of the highest fragment
		float fragmentBinSize=2.0f*(float)params.getFragmentTolerance().getTolerance(massPlusOne);
		double offset;
		
		if (fragmentBinSize>0.5f) {
			fragmentBinSize=lowResFragmentBinSize; // if tolerance is >0.25 Da, then jump to 1 Da to make use of the average amino acid mass defect
			offset=lowResFragmentBinOffset;
		} else if (fragmentBinSize<0.01f) {
			fragmentBinSize=0.01f;
			offset=0.0;
		} else {
			offset=0.0;
		}

		float inverseBinWidth=1.0f/fragmentBinSize;
		int arraySize=(int)((massPlusOne+fragmentBinSize+2.0)*inverseBinWidth);
		
		float[] binnedIntensityArray=new float[arraySize];
		int arraySizeMinusOne=arraySize-1;
		for (Peak peak : allPeaks) {
			int massIndex=(int)((peak.mass-offset)*inverseBinWidth);
			
			if (massIndex<0) massIndex=0;
			if (massIndex>=arraySize) massIndex=arraySize-1;
			if (binnedIntensityArray[massIndex]<peak.intensity) {
				binnedIntensityArray[massIndex]=peak.intensity;
			}
			
			// don't do this for low res fragment ions bin boundaries aren't an issue with the 0.4 offset
			if (fragmentBinSize<=0.5f&&addIntensityToNeighboringBins) {
				// neighboring intensities are 25 for b/y or 10 (the same) for neutral losses
				float neighboringIntensity=peak.intensity>neutralLossIntensity?peak.intensity/2.0f:peak.intensity;
				if (massIndex>0) {
					binnedIntensityArray[massIndex-1]=neighboringIntensity;
				}
				if (massIndex<arraySizeMinusOne) {
					binnedIntensityArray[massIndex+1]=neighboringIntensity;
				}
			}
		}
		return binnedIntensityArray;
	}
}

