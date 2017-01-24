package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class XCorrCalculator {
	 // set 50 to be the maximum value, see pp 982 bottom right
	private static float maxIntensity=50f;
	
	// divide spectrum into 10 equal regions, see pp 982 bottom right
	private static int groups=10; 
	
	// remove 10-u window around precursor, see pp 979 mid left
	private static double precursorRemovalMargin=5.0;
	
	public static Spectrum getTheoreticalSpectrum(String modifiedSequence, byte charge, SearchParameters params) {
		FragmentationType type=params.getFragType();
		AminoAcidConstants aaConstants=params.getAAConstants();
		FragmentationModel model=new FragmentationModel(modifiedSequence, aaConstants);
		
		// enforces uniqueness
		ArrayList<Peak> allPeaks=new ArrayList<Peak>();
		switch (type) {
			case YONLY:
				FragmentIon[] yIons=model.getYIons();
				allPeaks.addAll(getPeaks(yIons, 0.0, 50.0f));
				allPeaks.addAll(getPeaks(yIons, MassConstants.neutronMass, 25.0f));
				allPeaks.addAll(getPeaks(yIons, -MassConstants.neutronMass, 25.0f));
				allPeaks.addAll(getPeaks(yIons, -MassConstants.nh3, 10.0f));
				allPeaks.addAll(getPeaks(yIons, -MassConstants.oh2, 10.0f));
				allPeaks.addAll(getPeaks(yIons, -MassConstants.co, 10.0f));
				allPeaks.addAll(getPeaks(yIons, -MassConstants.nh3+MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(yIons, -MassConstants.oh2+MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(yIons, -MassConstants.co+MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(yIons, -MassConstants.nh3-MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(yIons, -MassConstants.oh2-MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(yIons, -MassConstants.co-MassConstants.neutronMass, 10.0f));
				break;
				
			case CID:
				FragmentIon[] yIonsCID=model.getYIons();
				allPeaks.addAll(getPeaks(yIonsCID, 0.0, 50.0f));
				allPeaks.addAll(getPeaks(yIonsCID, MassConstants.neutronMass, 25.0f));
				allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.neutronMass, 25.0f));
				allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.nh3, 10.0f));
				allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.oh2, 10.0f));
				allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.co, 10.0f));
				allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.nh3+MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.oh2+MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.co+MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.nh3-MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.oh2-MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(yIonsCID, -MassConstants.co-MassConstants.neutronMass, 10.0f));
				
				FragmentIon[] bIonsCID=model.getBIons();
				allPeaks.addAll(getPeaks(bIonsCID, 0.0, 50.0f));
				allPeaks.addAll(getPeaks(bIonsCID, MassConstants.neutronMass, 25.0f));
				allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.neutronMass, 25.0f));
				allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.nh3, 10.0f));
				allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.oh2, 10.0f));
				allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.co, 10.0f));
				allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.nh3+MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.oh2+MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.co+MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.nh3-MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.oh2-MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(bIonsCID, -MassConstants.co-MassConstants.neutronMass, 10.0f));
				break;
				
			case ETD:
				FragmentIon[] cIonsCID=model.getCIons();
				allPeaks.addAll(getPeaks(cIonsCID, 0.0, 50.0f));
				allPeaks.addAll(getPeaks(cIonsCID, MassConstants.neutronMass, 25.0f));
				allPeaks.addAll(getPeaks(cIonsCID, -MassConstants.neutronMass, 25.0f));
				allPeaks.addAll(getPeaks(cIonsCID, -MassConstants.nh3, 10.0f));
				allPeaks.addAll(getPeaks(cIonsCID, -MassConstants.oh2, 10.0f));
				allPeaks.addAll(getPeaks(cIonsCID, -MassConstants.nh3+MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(cIonsCID, -MassConstants.oh2+MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(cIonsCID, -MassConstants.nh3-MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(cIonsCID, -MassConstants.oh2-MassConstants.neutronMass, 10.0f));
				
				FragmentIon[] zIonsCID=model.getCIons();
				allPeaks.addAll(getPeaks(zIonsCID, 0.0, 50.0f));
				allPeaks.addAll(getPeaks(zIonsCID, MassConstants.neutronMass, 50.0f)); // z+1
				allPeaks.addAll(getPeaks(zIonsCID, 2*MassConstants.neutronMass, 25.0f));
				allPeaks.addAll(getPeaks(zIonsCID, -MassConstants.neutronMass, 25.0f));
				allPeaks.addAll(getPeaks(zIonsCID, -MassConstants.nh3, 10.0f));
				allPeaks.addAll(getPeaks(zIonsCID, -MassConstants.oh2, 10.0f));
				allPeaks.addAll(getPeaks(zIonsCID, -MassConstants.nh3+MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(zIonsCID, -MassConstants.oh2+MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(zIonsCID, -MassConstants.nh3-MassConstants.neutronMass, 10.0f));
				allPeaks.addAll(getPeaks(zIonsCID, -MassConstants.oh2-MassConstants.neutronMass, 10.0f));
				break;
				
			default:
				throw new EncyclopediaException("Unknown fragmentation type ["+type+"]");
		}
		
		// enforces exact mass uniqueness (but not nearby mass uniqueness)
		Collections.sort(allPeaks);
		int index=1;
		while (index<allPeaks.size()) {
			Peak current=allPeaks.get(index-1);
			Peak nextPeak=allPeaks.get(index);
			if (current.mass==nextPeak.mass) {
				if (nextPeak.intensity>current.intensity) {
					allPeaks.remove(index-1);
				} else {
					allPeaks.remove(index);
				}
			} else {
				index++;
			}
		}
		Pair<double[], float[]> peakArrays=Peak.toArrays(allPeaks);
		
		return getNormalizedSpectrum(model, charge, peakArrays.x, peakArrays.y);
	}
	
	public static ArrayList<Peak> getPeaks(FragmentIon[] ions, double delta, float intensity) {
		ArrayList<Peak> peaks=new ArrayList<Peak>();
		for (int i=0; i<ions.length; i++) {
			peaks.add(new Peak(ions[i].mass+delta, intensity));
		}
		return peaks;
	}
	
	
	/**
	 * see Eng et al, JASMS 1994
	 * @param s
	 * @param precursorMz
	 * @return
	 */
	public static Spectrum normalize(Spectrum s, double precursorMz) {
		double[] masses=s.getMassArray();
		float[] intensities=s.getIntensityArray();
		if (masses.length==0)
			return s;
		if (masses.length==1) {
			return getNormalizedSpectrum(s, masses, new float[] {maxIntensity});
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
		
		TDoubleArrayList normalizedMasses=new TDoubleArrayList();
		TFloatArrayList normalizedIntensities=new TFloatArrayList();
		
		General.divide(binMaxIntensity, maxIntensity);
		
		currentIndex=0;
		for (int i=0; i<intensities.length; i++) {
			if (masses[i]>minimumPrecursorRemoved&&masses[i]<maximumPrecursorRemoved) {
				continue;
			}
			
			while (masses[i]>binMaxMass[currentIndex]) {
				currentIndex++;
			}
			
			normalizedMasses.add(masses[i]);
			normalizedIntensities.add(intensities[i]/binMaxIntensity[currentIndex]);
		}
		
		return getNormalizedSpectrum(s, normalizedMasses.toArray(), normalizedIntensities.toArray());
	}

	public static Spectrum getNormalizedSpectrum(FragmentationModel model, byte precursorCharge, final double[] masses, final float[] intensities) {
		final float tic=General.sum(intensities);
		final String peptideModSeq=model.toString();
		final double mz=model.getChargedMass(precursorCharge);
		
		return getNormalizedSpectrum(masses, intensities, tic, 0.0f, peptideModSeq, mz);
	}
	
	public static Spectrum getNormalizedSpectrum(final Spectrum s, final double[] masses, final float[] intensities) {
		float tic=s.getTIC();
		float scanStartTime=s.getScanStartTime();
		double mz=s.getPrecursorMZ();
		String name=s.getSpectrumName();
		return getNormalizedSpectrum(masses, intensities, tic, scanStartTime, name, mz);
	}

	static Spectrum getNormalizedSpectrum(final double[] masses, final float[] intensities, final float tic, final float scanStartTime, final String name, final double mz) {
		return new Spectrum() {
			@Override
			public float getTIC() {
				return tic;
			}

			@Override
			public String getSpectrumName() {
				return name;
			}

			@Override
			public float getScanStartTime() {
				return scanStartTime;
			}

			@Override
			public double getPrecursorMZ() {
				return mz;
			}

			@Override
			public double[] getMassArray() {
				return masses;
			}

			@Override
			public float[] getIntensityArray() {
				return intensities;
			}
		};
	}
}
