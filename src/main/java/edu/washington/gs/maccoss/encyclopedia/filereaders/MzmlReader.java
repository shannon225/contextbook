package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import uk.ac.ebi.jmzml.model.mzml.BinaryDataArrayList;
import uk.ac.ebi.jmzml.model.mzml.CVParam;
import uk.ac.ebi.jmzml.model.mzml.Precursor;
import uk.ac.ebi.jmzml.model.mzml.PrecursorList;
import uk.ac.ebi.jmzml.model.mzml.Spectrum;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class MzmlReader {
	public static void main(String[] args) throws IOException, SQLException {
		StripeFile stripeFile=new StripeFile();
		stripeFile.openFile();

		Logger.logLine("Starting ...");
		File xmlFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/82593_lv_mcx_DIA_5mz_400to525.mzML");
		File saveFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/82593_lv_mcx_DIA_5mz_400to525.dia");
		
		MzMLUnmarshaller unmarshaller=new MzMLUnmarshaller(xmlFile);

		stripeFile.setFileName(xmlFile.getName(), unmarshaller.getMzMLId(), xmlFile.getAbsolutePath());

		int spectrumCount=unmarshaller.getObjectCountForXpath("/run/spectrumList/spectrum");
		Logger.logLine("Number of spectrum elements: "+spectrumCount);

		MzMLObjectIterator<Spectrum> spectrumIterator=unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", Spectrum.class);

		ArrayList<PrecursorScan> precursors=new ArrayList<PrecursorScan>();
		ArrayList<Stripe> stripes=new ArrayList<Stripe>();
		int count=0;
		int previousReport=0;
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

			String spectrumName=spectrum.getId();
			int spectrumIndex=spectrum.getIndex();
			HashMap<String, CVParam> cvparams=asCVMap(spectrum.getScanList().getScan().get(0).getCvParam());
			float scanStartTime=Float.parseFloat(cvparams.get("MS:1000016").getValue());
			BinaryDataArrayList bdal=spectrum.getBinaryDataArrayList();

			double[] massArray=ByteConverter.toDoubleArray(bdal.getBinaryDataArray().get(0).getBinaryDataAsNumberArray());
			float[] intensityArray=ByteConverter.toFloatArray(bdal.getBinaryDataArray().get(1).getBinaryDataAsNumberArray());

			if (p==null) {
				precursors.add(new PrecursorScan(spectrumName, spectrumIndex, scanStartTime, massArray, intensityArray));
			} else {
				HashMap<String, CVParam> isolationCVParams=asCVMap(p.getIsolationWindow().getCvParam());
				float isolationWindowTarget=Float.parseFloat(isolationCVParams.get("MS:1000827").getValue());
				float isolationWindowLowerOffset=Float.parseFloat(isolationCVParams.get("MS:1000828").getValue());
				float isolationWindowUpperOffset=Float.parseFloat(isolationCVParams.get("MS:1000829").getValue());
				
				stripes.add(new Stripe(spectrumName, p.getSpectrumRef(), spectrumIndex, scanStartTime, isolationWindowTarget-isolationWindowLowerOffset, isolationWindowTarget+isolationWindowUpperOffset, massArray, intensityArray));
			}
			
			if (precursors.size()>100) {
				stripeFile.addPrecursor(precursors);
				precursors.clear();
			}
			
			if (stripes.size()>100) {
				stripeFile.addStripe(stripes);
				stripes.clear();
			}
			int percent=(100*count)/spectrumCount;
			if (percent>previousReport) {
				previousReport=percent;
				Logger.logLine(percent+"% complete");
			}
			count++;
		}
		stripeFile.addPrecursor(precursors);
		stripeFile.addStripe(stripes);
		
		stripeFile.saveAsFile(saveFile);
		Logger.logLine("... Finished!");
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
