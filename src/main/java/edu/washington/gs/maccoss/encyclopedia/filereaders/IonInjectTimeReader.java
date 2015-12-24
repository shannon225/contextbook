package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import uk.ac.ebi.jmzml.model.mzml.CVParam;
import uk.ac.ebi.jmzml.model.mzml.Precursor;
import uk.ac.ebi.jmzml.model.mzml.PrecursorList;
import uk.ac.ebi.jmzml.model.mzml.Spectrum;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;

public class IonInjectTimeReader {
	public static void main(String[] args) throws IOException, SQLException {

		Logger.logLine("Starting ...");
		//File xmlFile=new File("/Users/searleb/Documents/freezer_experiment/freezer_data/121115_bcs_hela_24mz_400_1000_m20c_4D_151221163204_potentially_bad.mzML");
		File xmlFile=new File("/Users/searleb/Documents/freezer_experiment/freezer_data/121115_BCS_HeLa_24mz_400_1000.mzML");
		
		MzMLUnmarshaller unmarshaller=new MzMLUnmarshaller(xmlFile);

		int spectrumCount=unmarshaller.getObjectCountForXpath("/run/spectrumList/spectrum");
		Logger.logLine("Number of spectrum elements: "+spectrumCount);

		MzMLObjectIterator<Spectrum> spectrumIterator=unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", Spectrum.class);

		int count=0;
		int previousReport=0;
		
		HashMap<Range, ArrayList<XYPoint>> tracesByStripe=new HashMap<Range, ArrayList<XYPoint>>();
		while (spectrumIterator.hasNext()) {
			Spectrum spectrum=spectrumIterator.next();
			PrecursorList pl=spectrum.getPrecursorList();
			Precursor p=null;
			if (pl!=null) {
				for (Precursor precursor : pl.getPrecursor()) {
					p=precursor;
					break;
				}
			}

			HashMap<String, CVParam> cvparams=asCVMap(spectrum.getScanList().getScan().get(0).getCvParam());
			CVParam scanStartTimeCVParams=cvparams.get("MS:1000016");
			float multiplier;
			String unit=scanStartTimeCVParams.getUnitName();
			if ("second".equalsIgnoreCase(unit)) {
				multiplier=1.0f;
			} else if ("minute".equalsIgnoreCase(unit)) {
				multiplier=60.0f;
			} else if ("hour".equalsIgnoreCase(unit)) {
				multiplier=360.0f;
			} else if ("millisecond".equalsIgnoreCase(unit)) {
				multiplier=0.001f;
			} else {
				throw new EncyclopediaException("Unexpected time unit: "+unit);
			}
			
			float scanStartTime=multiplier*Float.parseFloat(scanStartTimeCVParams.getValue());

			
			CVParam ionInjectTimeCVParams=cvparams.get("MS:1000927");
			float ionInjectTime=Float.parseFloat(ionInjectTimeCVParams.getValue());

			Range range=null;
			if (p==null) {
				range=new Range(0, 2000);
			} else {
				HashMap<String, CVParam> isolationCVParams=asCVMap(p.getIsolationWindow().getCvParam());
				float isolationWindowTarget=Float.parseFloat(isolationCVParams.get("MS:1000827").getValue());
				float isolationWindowLowerOffset=Float.parseFloat(isolationCVParams.get("MS:1000828").getValue());
				float isolationWindowUpperOffset=Float.parseFloat(isolationCVParams.get("MS:1000829").getValue());
				
				//range=new Range(isolationWindowTarget-isolationWindowLowerOffset, isolationWindowTarget+isolationWindowUpperOffset);
			}

			if (range!=null) {
				ArrayList<XYPoint> stripe=tracesByStripe.get(range);
				if (stripe==null) {
					stripe=new ArrayList<XYPoint>();
					tracesByStripe.put(range, stripe);
				}
				stripe.add(new XYPoint(scanStartTime/60f, ionInjectTime));
			}
			
			int percent=(100*count)/spectrumCount;
			if (percent>previousReport) {
				previousReport=percent;
				Logger.logLine(percent+"% complete");
			}
			count++;
		}
		Logger.logLine("... Finished!");
		
		ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
		for (Entry<Range, ArrayList<XYPoint>> entry : tracesByStripe.entrySet()) {
			traces.add(new XYTrace(entry.getValue(), GraphType.line, entry.getKey().toString()));
		}
		
		Charter.launchChart("RT (min)", "Ion Inject Time", true, traces.toArray(new XYTrace[traces.size()]));
	}

	public static HashMap<String, CVParam> asCVMap(List<CVParam> params) {
		HashMap<String, CVParam> map=new HashMap<String, CVParam>();
		if (params==null) return map;
		for (CVParam cvParam : params) {
			map.put(cvParam.getAccession(), cvParam);
		}
		return map;
	}
}
