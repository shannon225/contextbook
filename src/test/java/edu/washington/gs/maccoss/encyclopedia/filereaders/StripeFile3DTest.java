package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYZPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYZTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.map.hash.TFloatObjectHashMap;
import gnu.trove.procedure.TFloatFloatProcedure;
import gnu.trove.procedure.TFloatObjectProcedure;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class StripeFile3DTest {
	public static void main(String[] args) throws Exception {
		File diaFile=new File("/Users/searleb/Documents/villen_manuscript/q05498_bs_MCF7_IMAC_DIA_R1.mzML");
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, PecanParameterParser.getDefaultParametersObject());
		
		List<PrecursorScan> scans=stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE);
		float rtResolution=10f;
		float mzResolution=10f;
		
		final TFloatObjectHashMap<TFloatFloatHashMap> threeDMap=new TFloatObjectHashMap<TFloatFloatHashMap>();
		
		for (PrecursorScan scan : scans) {
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
							points.add(new XYZPoint(rt, arg0, arg1));
							return true;
						}
					});
				}
				return true;
			}
		});
		
		System.out.println(points.size());

		Charter.launchChart("Retention Time (Seconds)", "M/z", false, new XYZTrace("Intensity", points));
	}

}
