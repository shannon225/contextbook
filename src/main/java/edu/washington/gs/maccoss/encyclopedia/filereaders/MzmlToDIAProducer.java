package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.set.hash.TIntHashSet;
import uk.ac.ebi.jmzml.model.mzml.BinaryDataArrayList;
import uk.ac.ebi.jmzml.model.mzml.CVParam;
import uk.ac.ebi.jmzml.model.mzml.Precursor;
import uk.ac.ebi.jmzml.model.mzml.PrecursorList;
import uk.ac.ebi.jmzml.model.mzml.Spectrum;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;

public class MzmlToDIAProducer implements Runnable {
	private final BlockingQueue<MzmlBlock> mzmlBlockQueue;
	private final MzMLUnmarshaller unmarshaller;
	private final SearchParameters parameters;
	private final HashMap<Range, TFloatArrayList> retentionTimesByStripe=new HashMap<Range, TFloatArrayList>();

	public MzmlToDIAProducer(MzMLUnmarshaller unmarshaller, BlockingQueue<MzmlBlock> mzmlBlockQueue, SearchParameters parameters) {
		this.unmarshaller=unmarshaller;
		this.mzmlBlockQueue=mzmlBlockQueue;
		this.parameters=parameters;
	}
	
	public HashMap<Range, TFloatArrayList> getRetentionTimesByStripe() {
		return retentionTimesByStripe;
	}

	@Override
	public void run() {
		try {
		int spectrumCount=unmarshaller.getObjectCountForXpath("/run/spectrumList/spectrum");
		Logger.logLine("Number of spectrum elements: "+spectrumCount);

		MzMLObjectIterator<Spectrum> spectrumIterator=unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", Spectrum.class);

		ArrayList<PrecursorScan> precursors=new ArrayList<PrecursorScan>();
		ArrayList<Stripe> stripes=new ArrayList<Stripe>();
		int count=0;
		int previousReport=0;
		
		float defaultOffset=parameters.getPrecursorWindowSize()/2.0f;
		TIntHashSet spectrumIndices=new TIntHashSet();
		
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
			if (spectrumIndices.contains(spectrumIndex)) {
				System.out.println("Duplicate\t"+spectrumIndex);
			} else {
				spectrumIndices.add(spectrumIndex);
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
			BinaryDataArrayList bdal=spectrum.getBinaryDataArrayList();
			if (bdal==null||bdal.getBinaryDataArray()==null||bdal.getBinaryDataArray().get(0)==null||bdal.getBinaryDataArray().get(1)==null) {
				Logger.errorLine("Skipping unexpected empty binary data for spectrum"+" ("+spectrumIndex+"): "+spectrumName);
				continue;
			}

			double[] massArray=ByteConverter.toDoubleArray(bdal.getBinaryDataArray().get(0).getBinaryDataAsNumberArray());
			float[] intensityArray=ByteConverter.toFloatArray(bdal.getBinaryDataArray().get(1).getBinaryDataAsNumberArray());

			if (p==null) {
				if (parameters.getPrecursorOffsetPPM()!=0.0) {
					double[] deltaArray=General.multiply(massArray, parameters.getPrecursorOffsetPPM()/1000000.0);
					massArray=General.subtract(massArray, deltaArray);
				}
				precursors.add(new PrecursorScan(spectrumName, spectrumIndex, scanStartTime, massArray, intensityArray));
			} else {
				HashMap<String, CVParam> isolationCVParams=asCVMap(p.getIsolationWindow().getCvParam());
				float isolationWindowTarget=Float.parseFloat(isolationCVParams.get("MS:1000827").getValue());
				CVParam lowerParam=isolationCVParams.get("MS:1000828");
				CVParam upperParam=isolationCVParams.get("MS:1000829");
				float isolationWindowLowerOffset;
				float isolationWindowUpperOffset;
				if (lowerParam==null||upperParam==null) {
					if (defaultOffset<=0) {
						throw new EncyclopediaException("Error reading mzML! Precursor window offsets not specified and no default window size specified!");
					}
					isolationWindowLowerOffset=defaultOffset;
					isolationWindowUpperOffset=defaultOffset;
				} else {
					isolationWindowLowerOffset=Float.parseFloat(lowerParam.getValue());
					isolationWindowUpperOffset=Float.parseFloat(upperParam.getValue());
				}

				if (parameters.getFragmentOffsetPPM()!=0.0) {
					double[] deltaArray=General.multiply(massArray, parameters.getFragmentOffsetPPM()/1000000.0);
					massArray=General.subtract(massArray, deltaArray);
				}
				Stripe stripe=new Stripe(spectrumName, p.getSpectrumRef(), spectrumIndex, scanStartTime, isolationWindowTarget-isolationWindowLowerOffset, isolationWindowTarget+isolationWindowUpperOffset, massArray, intensityArray);
				stripes.add(stripe);
				Range range=stripe.getRange();
				TFloatArrayList stripeRTs=retentionTimesByStripe.get(range);
				if (stripeRTs==null) {
					stripeRTs=new TFloatArrayList();
					retentionTimesByStripe.put(range, stripeRTs);
				}
				stripeRTs.add(scanStartTime);
			}

			try {
				if (precursors.size()>100||stripes.size()>1000) {
					mzmlBlockQueue.put(new MzmlBlock(precursors, stripes));
					precursors.clear();
					stripes.clear();
				}
			} catch (InterruptedException ie) {
				Logger.errorLine("Mzml reading interrupted!");
				Logger.errorException(ie);
			}
			int percent=(100*count)/spectrumCount;
			if (percent>previousReport) {
				previousReport=percent;
				Logger.logLine(percent+"% complete");
			}
			count++;
		}
		try {
			mzmlBlockQueue.put(new MzmlBlock(precursors, stripes));
			mzmlBlockQueue.put(MzmlBlock.POISON_BLOCK);			
		} catch (InterruptedException ie) {
			Logger.errorLine("Mzml reading interrupted!");
			Logger.errorException(ie);
		}
		} catch (Exception e) {
			Logger.errorLine("Mzml reading failed!");
			Logger.errorException(e);
		}
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
