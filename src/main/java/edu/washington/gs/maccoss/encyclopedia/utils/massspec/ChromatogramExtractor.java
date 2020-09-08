package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;

public class ChromatogramExtractor {
	public static final byte[] isotopes=new byte[] {0, 1, 2};
	public static final Color[] isotopeColors=new Color[] {new Color(0, 0, 255), new Color(138, 43, 226), new Color(165, 42, 42)};
	public static XYTraceInterface[] extractPrecursorChromatograms(MassTolerance tolerance, double precursorMz, byte charge, List<? extends Spectrum> precursors) {
		double[] targetMasses=new double[isotopes.length];
		for (int i=0; i<targetMasses.length; i++) {
			targetMasses[i] = MassConstants.getChargedIsotopeMass(precursorMz, charge, isotopes[i]);
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
			XYTraceInterface trace=new XYTrace(traces[i], GraphType.line, name, isotopeColors[i], 3.0f);
			kept.add(SkylineSGFilter.paddedSavitzkyGolaySmooth(trace));
		}
		return kept.toArray(new XYTrace[kept.size()]);
	}

	public static HashMap<FragmentIon, XYTrace> extractFragmentChromatograms(MassTolerance tolerance, FragmentIon[] ionTypes, List<? extends Spectrum> stripes, Float targetRTInSec, GraphType type) {
		HashMap<FragmentIon, XYTrace> kept=new HashMap<FragmentIon, XYTrace>();

		ArrayList<FragmentIon> centerIonTypes=new ArrayList<FragmentIon>();
		
		if (targetRTInSec==null) {
			centerIonTypes.addAll(Arrays.asList(ionTypes));
		} else {
			Spectrum bestStripe=getTargetStripeByRT(stripes, targetRTInSec);
			// no signal of any kind at retention time!
			if (bestStripe==null) return kept;

			for (int i=0; i<ionTypes.length; i++) {
				float intensity=tolerance.getIntegratedIntensity(bestStripe.getMassArray(), bestStripe.getIntensityArray(), ionTypes[i].getMass());
				if (intensity>0) {
					centerIonTypes.add(ionTypes[i]);
				}
			}
		}
		
		HashMap<FragmentIon, ArrayList<XYPoint>> traces=new HashMap<FragmentIon, ArrayList<XYPoint>>();
		for (FragmentIon centerIon : centerIonTypes) {
			traces.put(centerIon, new ArrayList<XYPoint>());
		}

		for (int i=0; i<stripes.size(); i++) {
			Spectrum spectrum=stripes.get(i);
			double[] massArray=spectrum.getMassArray();
			float[] intensityArray=spectrum.getIntensityArray();
			for (FragmentIon centerIon : centerIonTypes) {
				float intensity=tolerance.getIntegratedIntensity(massArray, intensityArray, centerIon.getMass());
				ArrayList<XYPoint> points=traces.get(centerIon);
				points.add(new XYPoint(spectrum.getScanStartTime()/60, intensity));
			}
		}
		
		for (Entry<FragmentIon, ArrayList<XYPoint>> traceData : traces.entrySet()) {
			FragmentIon key=traceData.getKey();
			String name=key.toString();
			XYTrace trace=null;
			switch (type) {
			case boldline:
				trace=new XYTrace(traceData.getValue(), GraphType.line, name, key.getColor(), 3.0f);
				break;
			case line:
				trace=new XYTrace(traceData.getValue(), GraphType.line, name, key.getColor(), 2.0f);
				break;
			case dashedline:
				trace=new XYTrace(traceData.getValue(), GraphType.dashedline, name, key.getColor(), 1.0f);
				break;
			default:
				trace=new XYTrace(traceData.getValue(), GraphType.line, name, key.getColor(), 2.0f);
				break;
			}
			XYTrace sgSmoothed=SkylineSGFilter.paddedSavitzkyGolaySmooth(trace);
			kept.put(key, sgSmoothed);
			//kept.put(key, trace);
		}
		
		return kept;
	}

	public static Spectrum getTargetStripeByRT(List<? extends Spectrum> stripes, Float targetRTInSec) {
		Spectrum bestStripe=null;
		float bestDelta=Float.MAX_VALUE;
		for (Spectrum stripe : stripes) {
			float delta=Math.abs(stripe.getScanStartTime()-targetRTInSec);
			if (delta<bestDelta) {
				bestDelta=delta;
				bestStripe=stripe;
			}
		}
		return bestStripe;
	}

}
