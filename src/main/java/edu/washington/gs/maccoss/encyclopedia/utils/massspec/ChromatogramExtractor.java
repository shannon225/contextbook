package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;

public class ChromatogramExtractor {
	public static final byte[] isotopes=new byte[] {0, 1, 2};
	public static final Color[] isotopeColors=new Color[] {Color.BLUE, Color.RED, new Color(0, 180, 0)};
	public static XYTrace[] extractPrecursorChromatograms(MassTolerance tolerance, double precursorMz, byte charge, ArrayList<Spectrum> precursors) {
		double[] targetMasses=new double[isotopes.length];
		for (int i=0; i<targetMasses.length; i++) {
			targetMasses[i]=precursorMz+(isotopes[i]*MassConstants.neutronMass/charge);
		}
		@SuppressWarnings("unchecked")
		ArrayList<XYPoint>[] traces=new ArrayList[targetMasses.length];
		for (int i=0; i<traces.length; i++) {
			traces[i]=new ArrayList<XYPoint>();
		}
		
		for (Spectrum spectrum : precursors) {
			double[] massArray=spectrum.getMassArray();
			float[] intensityArray=spectrum.getIntensityArray();
			
			for (int i=0; i<targetMasses.length; i++) {
				float intensity=tolerance.getIntegratedIntensity(massArray, intensityArray, targetMasses[i]);
				traces[i].add(new XYPoint(spectrum.getScanStartTime()/60, intensity));
			}
		}
		ArrayList<XYTrace> kept=new ArrayList<XYTrace>();
		for (int i=0; i<traces.length; i++) {
			String name;
			if (isotopes[i]>0) {
				name="Precursor+"+isotopes[i];
			} else {
				name="Precursor";
			}
			XYTrace trace=new XYTrace(traces[i], GraphType.line, name, isotopeColors[i], 3.0f);
			kept.add(SkylineSGFilter.paddedSavitzkyGolaySmooth(trace));
		}
		return kept.toArray(new XYTrace[kept.size()]);
	}
	
	public static HashMap<FragmentIon, XYTrace> extractFragmentChromatograms(MassTolerance tolerance, FragmentIon[] ionTypes, ArrayList<Spectrum> stripes, float targetRTInSec) {
		@SuppressWarnings("unchecked")
		ArrayList<XYPoint>[] traces=new ArrayList[ionTypes.length];
		boolean[] gotTrace=new boolean[traces.length];
		for (int i=0; i<traces.length; i++) {
			traces[i]=new ArrayList<XYPoint>();
		}
		for (int i=0; i<stripes.size(); i++) {
			Spectrum spectrum=stripes.get(i);
			double[] massArray=spectrum.getMassArray();
			float[] intensityArray=spectrum.getIntensityArray();
			for (int j=0; j<ionTypes.length; j++) {
				if (ionTypes[j].index<=2) continue;
				
				float intensity=tolerance.getIntegratedIntensity(massArray, intensityArray, ionTypes[j].mass);
				traces[j].add(new XYPoint(spectrum.getScanStartTime()/60, intensity));
				if (intensity>0) gotTrace[j]=true;
			}
		}
		
		// FIXME REMOVE IONS NOT IN CENTER!
		
		HashMap<FragmentIon, XYTrace> kept=new HashMap<FragmentIon, XYTrace>();
		for (int i=0; i<traces.length; i++) {
			if (gotTrace[i]) {
				String name=ionTypes[i].toString();
				XYTrace trace=new XYTrace(traces[i], GraphType.line, name, RandomGenerator.randomColor(name.hashCode()), 3.0f);
				kept.put(ionTypes[i], SkylineSGFilter.paddedSavitzkyGolaySmooth(trace));
			}
		}
		return kept;
	}

}
