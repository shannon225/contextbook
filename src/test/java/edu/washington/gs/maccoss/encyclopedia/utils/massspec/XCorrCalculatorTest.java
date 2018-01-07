package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.awt.Dimension;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import junit.framework.TestCase;

public class XCorrCalculatorTest extends TestCase {
	private static final SearchParameters MAIN_PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.HCD, new MassTolerance(10, MassErrorUnitType.PPM), new MassTolerance(10, MassErrorUnitType.PPM), DigestionEnzyme.getEnzyme("trypsin"), false, true);
	//private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(0.05, MassErrorUnitType.AMU), new MassTolerance(10, MassErrorUnitType.PPM), DigestionEnzyme.getEnzyme("trypsin"));
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(1, MassErrorUnitType.AMU), new MassTolerance(1, MassErrorUnitType.AMU), DigestionEnzyme.getEnzyme("trypsin"), false, true);
	//private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.YONLY, new MassTolerance(10, MassErrorUnitType.PPM), new MassTolerance(10, MassErrorUnitType.PPM), DigestionEnzyme.getEnzyme("trypsin"));
	
	public static void main3(String[] args) {
		// timing test
		final byte charge=2;
		final float chargedMz=(float)((1329.6335+(charge-1)*MassConstants.protonMass)/charge);
		
		Spectrum s=getSDFHLFGPPGKK();

		SparseXCorrCalculator preprocessedmodel=new SparseXCorrCalculator("SDFHLFGPPGKK", charge, MAIN_PARAMETERS);
		SparseXCorrSpectrum sparse=preprocessedmodel.normalize(s, new Range(chargedMz-10f, chargedMz+10f));

		long time=System.currentTimeMillis();
		for (int i=0; i<100000; i++) {
			preprocessedmodel.score(sparse);	
		}
		System.out.println("100000 Cached Sparse: "+(System.currentTimeMillis()-time));

		ArrayXCorrCalculator arraypreprocessedmodel=new ArrayXCorrCalculator("SDFHLFGPPGKK", chargedMz, charge, MAIN_PARAMETERS);
		float[] normalized=arraypreprocessedmodel.normalize(s);
		
		time=System.currentTimeMillis();
		for (int i=0; i<100000; i++) {
			arraypreprocessedmodel.score(normalized);
		}
		System.out.println("100000 Cached Array: "+(System.currentTimeMillis()-time));
		
		time=System.currentTimeMillis();
		for (int i=0; i<1000; i++) {
			preprocessedmodel=new SparseXCorrCalculator("SDFHLFGPPGKK", charge, MAIN_PARAMETERS);
			sparse=preprocessedmodel.normalize(s, new Range(chargedMz-10f, chargedMz+10f));
			preprocessedmodel.score(sparse);	
		}
		System.out.println("Sparse: "+(System.currentTimeMillis()-time));

		
		time=System.currentTimeMillis();
		for (int i=0; i<1000; i++) {
			arraypreprocessedmodel=new ArrayXCorrCalculator("SDFHLFGPPGKK", chargedMz, charge, MAIN_PARAMETERS);
			normalized=arraypreprocessedmodel.normalize(s);
			arraypreprocessedmodel.score(normalized);
		}
		System.out.println("Array: "+(System.currentTimeMillis()-time));
	}
	
	public static void main(String[] args) {
		final byte charge=2;
		final float chargedMz=(float)((1329.6335+(charge-1)*MassConstants.protonMass)/charge);
		
		Spectrum s=getSDFHLFGPPGKK();

		SparseXCorrSpectrum f=SparseXCorrCalculator.normalize(s, new Range(chargedMz-10.0f, chargedMz+10.0f), false, PARAMETERS);
		SparseXCorrSpectrum t=SparseXCorrCalculator.getTheoreticalSpectrum("SDFHLFGPPGKK", charge, PARAMETERS);
		Charter.launchChart(getNormalizedSpectrum(s, SparseXCorrCalculator.biggestFragmentMass, charge, f, PARAMETERS), "Spectrum", new Dimension(300, 200));
		Charter.launchChart(getNormalizedSpectrum(s, SparseXCorrCalculator.biggestFragmentMass, charge, t, PARAMETERS), "Model", new Dimension(300, 200));
		f=SparseXCorrCalculator.preprocessSpectrum(f);
		Charter.launchChart(getNormalizedSpectrum(s, SparseXCorrCalculator.biggestFragmentMass, charge, f, PARAMETERS), "PP Spectrum", new Dimension(300, 200));
		t=SparseXCorrCalculator.preprocessSpectrum(t);
		Charter.launchChart(getNormalizedSpectrum(s, SparseXCorrCalculator.biggestFragmentMass, charge, t, PARAMETERS), "PP Model", new Dimension(300, 200));
	}
	
	public void testXCorr() {
		final byte charge=2;
		final float chargedMz=(float)((1329.6335+(charge-1)*MassConstants.protonMass)/charge);
		
		Spectrum s=getSDFHLFGPPGKK();
		
		SparseXCorrCalculator preprocessedSpectrum=new SparseXCorrCalculator(s, new Range(chargedMz-10.0f, chargedMz+10.0f), PARAMETERS);
		float spectrumFirst=preprocessedSpectrum.score("SDFHLFGPPGKK", charge);
		System.out.println("spectrumFirst xcorr: "+spectrumFirst);
		
		SparseXCorrCalculator preprocessedmodel=new SparseXCorrCalculator("SDFHLFGPPGKK", charge, PARAMETERS);
		float modelFirst=preprocessedmodel.score(s, new Range(chargedMz-10.0f, chargedMz+10.0f));
		System.out.println("modelFirst xcorr: "+modelFirst);
		
		
		SparseXCorrSpectrum f=SparseXCorrCalculator.normalize(s, new Range(chargedMz-10.0f, chargedMz+10.0f), false, PARAMETERS);
		SparseXCorrSpectrum t=SparseXCorrCalculator.getTheoreticalSpectrum("SDFHLFGPPGKK", charge, PARAMETERS);

		for (int i=-75; i<=75; i++) {
			System.out.println(i+"\t"+t.dotProduct(f, i));
		}
		
		float center=t.dotProduct(f, 0);
		float avg=0.0f;
		for (int i=-75; i<75; i++) {
			if (i!=0) {
				avg+=t.dotProduct(f, i);
			}
		}
		avg=avg/150.0f;
		
		float originalCalculation=(center-avg)/1e4f;
		

		float center2=f.dotProduct(t, 0);
		float avg2=0.0f;
		for (int i=-75; i<75; i++) {
			if (i!=0) {
				avg2+=f.dotProduct(t, i);
			}
		}
		avg2=avg2/150.0f;
		
		float originalCalculation2=(center2-avg2)/1e4f;
		assertEquals(originalCalculation, originalCalculation2, 0.05f);
		assertEquals(originalCalculation, modelFirst, 0.05f);
		assertEquals(originalCalculation, spectrumFirst, 0.05f);
		System.out.println(center+"\t"+avg+"\t"+originalCalculation);

		ArrayXCorrCalculator arrayPreprocessedSpectrum=new ArrayXCorrCalculator(s, chargedMz, charge, PARAMETERS);
		float arraySpectrumFirst=arrayPreprocessedSpectrum.score("SDFHLFGPPGKK");
		System.out.println("array spectrumFirst xcorr: "+spectrumFirst);
		
		ArrayXCorrCalculator arrayPreprocessedmodel=new ArrayXCorrCalculator("SDFHLFGPPGKK", chargedMz, charge, PARAMETERS);
		float arrayModelFirst=arrayPreprocessedmodel.score(s);
		System.out.println("array modelFirst xcorr: "+modelFirst);
		assertEquals(originalCalculation, arraySpectrumFirst, 0.05f);
		assertEquals(originalCalculation, arrayModelFirst, 0.05f);
		
		//s=getNormalizedSpectrum(s, chargedMz, charge, f, PARAMETERS);
		//Charter.launchChart(s, "model orig:"+originalCalculation2+" spec orig:"+originalCalculation+" model:"+modelFirst+" spec:"+spectrumFirst);
	}

	public static Spectrum getSDFHLFGPPGKK() { 
		// from TRFE_CHICK

		final byte charge=2;
		final double chargedMz=(1329.6335+(charge-1)*MassConstants.protonMass)/charge;
		
		InputStream is=XCorrCalculatorTest.class.getResourceAsStream("/040203_XXX_X1_1_OT_5seq.02.00085.2.dta.txt");
		ArrayList<Peak> peaks=getData(is);
		Pair<double[], float[]> peakArrays=Peak.toArrays(peaks);
		return getSpectrum(peakArrays.x, General.protectedSqrt(peakArrays.y), General.sum(peakArrays.y), 0.0f, "SDFHLFGPPGKK", chargedMz); 
	}

	public static ArrayList<Peak> getData(InputStream is) {
		final ArrayList<Peak> rts=new ArrayList<Peak>();

		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				double mass=Double.parseDouble(row.get("mass"));
				float intensity=Float.parseFloat(row.get("intensity"));
				rts.add(new Peak(mass, intensity));
			}
		};

		BlockingQueue<Map<String, String>> blockingQueue=new LinkedBlockingQueue<Map<String, String>>();
		TableParserProducer producer=new TableParserProducer(blockingQueue, is, " ", 1);
		TableParserConsumer consumer=new TableParserConsumer(blockingQueue, muscle);

		Thread producerThread=new Thread(producer);
		Thread consumerThread=new Thread(consumer);
		producerThread.start();
		consumerThread.start();

		try {
			producerThread.join();
			consumerThread.join();
		} catch (InterruptedException ie) {
			Logger.errorLine("Percolator reading interrupted!");
			Logger.errorException(ie);
		}

		return rts;
	}
	
	static Spectrum getNormalizedSpectrum(final Spectrum s, final double precursorMz, final byte charge, final SparseXCorrSpectrum intensityBins, final SearchParameters params) {
		double massPlusOne=precursorMz*charge-(charge-1)*MassConstants.protonMass;
		// set tolerance to 2x the fragment tolerance of the highest fragment
		float fragmentBinSize=2.0f*(float) params.getFragmentTolerance().getTolerance(massPlusOne);
		double offset;
		if (fragmentBinSize>0.5f) {
			fragmentBinSize=ArrayXCorrCalculator.lowResFragmentBinSize; 
			offset=ArrayXCorrCalculator.lowResFragmentBinOffset;
		} else if (fragmentBinSize<0.01f) {
			fragmentBinSize=0.01f;
			offset=0.0;
		} else {
			offset=0.0;
		}
		
		float tic=s.getTIC();
		float scanStartTime=s.getScanStartTime();
		double mz=s.getPrecursorMZ();
		String name=s.getSpectrumName();
		
		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();
		int[] indices=intensityBins.getIndices();
		float[] indexedIntensities=intensityBins.getIntensityArray();
		for (int i=0; i<indices.length; i++) {
			double mass=(indices[i]*fragmentBinSize)+offset;
			masses.add(mass);
			intensities.add(indexedIntensities[i]);
		}
		
		float[] intensityArray=intensities.toArray();
		return getSpectrum(masses.toArray(), intensityArray, tic, scanStartTime, name, mz);
	}

	static Spectrum getSpectrum(final double[] masses, final float[] intensities, final float tic, final float scanStartTime, final String name, final double mz) {
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
