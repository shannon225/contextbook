package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.map.hash.TIntFloatHashMap;
import gnu.trove.procedure.TIntFloatProcedure;

public class StripeFileTest {
	public static void main(String[] args) throws Exception {
		File diaFile=new File("/Users/searleb/Downloads/Pool_600.dia");
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, PecanParameterParser.getDefaultParametersObject());
		
		ArrayList<PrecursorScan> scans=stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE);
		ArrayList<XYPoint> tic=new ArrayList<XYPoint>();
		TIntFloatHashMap mzSum=new TIntFloatHashMap();
		for (PrecursorScan scan : scans) {
			float rt=scan.getScanStartTime()/60f;
			double[] masses=scan.getMassArray();
			float[] intensities=scan.getIntensityArray();
			float intensity=General.sum(intensities);
			XYPoint point=new XYPoint(rt, intensity);
			tic.add(point);
			
			for (int i=0; i<intensities.length; i++) {
				int index=(int)Math.round(masses[i]);
				mzSum.adjustOrPutValue(index, intensities[i], intensities[i]);
			}
		}

		final ArrayList<XYPoint> mzs=new ArrayList<XYPoint>();
		mzSum.forEachEntry(new TIntFloatProcedure() {
			@Override
			public boolean execute(int arg0, float arg1) {
				mzs.add(new XYPoint(arg0, arg1));
				System.out.println(arg0+"\t"+arg1);
				return true;
			}
		});

		//Charter.launchChart("Retention Time (Minutes)", "Total Ion Chromatogram", false, new XYTrace(tic, GraphType.line, "TIC"));
		Charter.launchChart("M/z", "Total Ion Chromatogram", false, new XYTrace(mzs, GraphType.line, "Mz"));
	}

}
