package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class XCorrCalculator {

	// set 50 to be the maximum value, see pp 982 bottom right
	private static float primaryIonIntensity=50.0f;
	 private static final float neutralLossIntensity=10.0f;
	
	// divide spectrum into 10 equal regions, see pp 982 bottom right
	private static int groups=10; 
	
	// remove 10-u window around precursor, see pp 979 mid left
	private static double precursorRemovalMargin=5.0;
	
	/**
	 * see Eng et al, JASMS 1994
	 * @param s
	 * @param precursorMz
	 * @return
	 */
	public static float[] normalize(Spectrum s, double precursorMz, byte charge, boolean addIntensityToNeighboringBins, SearchParameters params) {
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
		
		General.divide(binMaxIntensity, primaryIonIntensity);
		
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
	
	public static float[] getTheoreticalSpectrum(String modifiedSequence, double precursorMz, byte charge, SearchParameters params) {
		double massPlusOne=precursorMz*charge-(charge-1)*MassConstants.protonMass;
		
		FragmentationType type=params.getFragType();
		AminoAcidConstants aaConstants=params.getAAConstants();
		FragmentationModel model=new FragmentationModel(modifiedSequence, aaConstants);
		
		ArrayList<Peak> allPeaks=new ArrayList<Peak>();
		switch (type) {
			case YONLY:
				FragmentIon[] yIons=model.getYIons();
				allPeaks.addAll(getPeaks(yIons, 0.0, primaryIonIntensity));
				allPeaks.addAll(getPeaks(yIons, -MassConstants.nh3, neutralLossIntensity));
				allPeaks.addAll(getPeaks(yIons, -MassConstants.oh2, neutralLossIntensity));
				break;
				
			case CID:
				FragmentIon[] yIonsCID=model.getYIons();
				allPeaks.addAll(getPeaks(yIonsCID, 0.0, primaryIonIntensity));
				allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.nh3, neutralLossIntensity));
				allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.oh2, neutralLossIntensity));
				
				FragmentIon[] bIonsCID=model.getBIons();
				allPeaks.addAll(getPeaks(bIonsCID, 0.0, primaryIonIntensity));
				allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.nh3, neutralLossIntensity));
				allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.oh2, neutralLossIntensity));
				allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.co, neutralLossIntensity));
				break;
				
			case ETD:
				FragmentIon[] cIonsCID=model.getCIons();
				allPeaks.addAll(getPeaks(cIonsCID, 0.0, primaryIonIntensity));
				allPeaks.addAll(getPeaks(cIonsCID, -MassConstants.nh3, neutralLossIntensity));
				allPeaks.addAll(getPeaks(cIonsCID, -MassConstants.oh2, neutralLossIntensity));
				
				FragmentIon[] zIonsCID=model.getCIons();
				allPeaks.addAll(getPeaks(zIonsCID, 0.0, primaryIonIntensity));
				allPeaks.addAll(getPeaks(zIonsCID, MassConstants.neutronMass, primaryIonIntensity)); // z+1
				allPeaks.addAll(getPeaks(zIonsCID, -MassConstants.nh3, neutralLossIntensity));
				allPeaks.addAll(getPeaks(zIonsCID, -MassConstants.oh2, neutralLossIntensity));
				break;
				
			default:
				throw new EncyclopediaException("Unknown fragmentation type ["+type+"]");
		}
		
		return getIntensityArray(params, allPeaks, massPlusOne, true);
	}
	
	private static ArrayList<Peak> getPeaks(FragmentIon[] ions, double delta, float intensity) {
		ArrayList<Peak> peaks=new ArrayList<Peak>();
		for (int i=0; i<ions.length; i++) {
			peaks.add(new Peak(ions[i].mass+delta, intensity));
		}
		return peaks;
	}

	private static float[] getIntensityArray(SearchParameters params, ArrayList<Peak> allPeaks, double massPlusOne, boolean addIntensityToNeighboringBins) {
		Collections.sort(allPeaks);
		
		// set tolerance to 2x the fragment tolerance of the highest fragment
		float fragmentBinSize=2.0f*(float)params.getFragmentTolerance().getTolerance(massPlusOne);
		if (fragmentBinSize<0.01f) fragmentBinSize=0.01f; 
		float inverseBinWidth=1.0f/fragmentBinSize;
		int arraySize=(int)((massPlusOne+fragmentBinSize+2.0)*inverseBinWidth);
		
		float[] binnedIntensityArray=new float[arraySize];
		int arraySizeMinusOne=arraySize-1;
		for (Peak peak : allPeaks) {
			int massIndex=(int)(peak.mass*inverseBinWidth);
			if (massIndex>=arraySize) massIndex=arraySize-1;
			if (binnedIntensityArray[massIndex]<peak.intensity) {
				binnedIntensityArray[massIndex]=peak.intensity;
			}
			
			if (addIntensityToNeighboringBins) {
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
