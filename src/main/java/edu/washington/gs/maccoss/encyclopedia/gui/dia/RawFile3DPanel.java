package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.util.ArrayList;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYZPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYZTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.AcquiredSpectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.map.hash.TFloatObjectHashMap;
import gnu.trove.procedure.TFloatFloatProcedure;
import gnu.trove.procedure.TFloatObjectProcedure;

public class RawFile3DPanel {
	public static ChartPanel getPrecursorsPanel(ArrayList<PrecursorScan> precursors) {
		ArrayList<AcquiredSpectrum> spectra=new ArrayList<>();
		spectra.addAll(precursors);
		return getPanelInternal(spectra);
	}
	public static ChartPanel getFragmentsPanel(ArrayList<FragmentScan> fragments) {
		ArrayList<AcquiredSpectrum> spectra=new ArrayList<>();
		spectra.addAll(fragments);
		return getPanelInternal(spectra);
	}
	private static ChartPanel getPanelInternal(ArrayList<AcquiredSpectrum> precursors) {

		float rtResolution=10f;
		float mzResolution=10f;
		
		final TFloatObjectHashMap<TFloatFloatHashMap> threeDMap=new TFloatObjectHashMap<TFloatFloatHashMap>();
		
		for (AcquiredSpectrum scan : precursors) {
			final float rt=Math.round(scan.getScanStartTime()/rtResolution)*rtResolution;
			double[] masses=scan.getMassArray();
			float[] intensities=scan.getIntensityArray();
			
			TFloatFloatHashMap mzSum=threeDMap.get(rt);
			if (mzSum==null) {
				mzSum=new TFloatFloatHashMap();
				threeDMap.put(rt, mzSum);
			}
			for (int i=0; i<intensities.length; i++) {
				float mzIndex=(int)Math.round(masses[i]/mzResolution)*mzResolution;
				float intensity=Log.log10(intensities[i]);
				mzSum.adjustOrPutValue(mzIndex, intensity, intensity);
			}
		}

		final ArrayList<XYZPoint> points=new ArrayList<XYZPoint>();
		threeDMap.forEachEntry(new TFloatObjectProcedure<TFloatFloatHashMap>() {
			@Override
			public boolean execute(final float rt, TFloatFloatHashMap arg1) {
				TFloatFloatHashMap mzSum=threeDMap.get(rt);
				if (mzSum!=null) {
					mzSum.forEachEntry(new TFloatFloatProcedure() {
						@Override
						public boolean execute(float arg0, float arg1) {
							points.add(new XYZPoint(rt/60f, arg0, arg1));
							return true;
						}
					});
				}
				return true;
			}
		});
		
		System.out.println(points.size());

		return Charter.getChart("Retention Time (min)", "m/z", false, new XYZTrace("Intensity", points));
	}

}
