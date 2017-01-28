package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.MedianInterpolatorTest;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import junit.framework.TestCase;

public class XCorrCalculatorTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(0.5, MassErrorUnitType.AMU), new MassTolerance(10, MassErrorUnitType.PPM), DigestionEnzyme.getEnzyme("trypsin"));

	public static void main(String[] args) {
		final byte charge=2;
		final double chargedMz=(1329.6335+(charge-1)*MassConstants.protonMass)/charge;
		
		Spectrum s=getSDFHLFGPPGKK();

		float[] f=XCorrCalculator.normalize(s, s.getPrecursorMZ(), charge, false, PARAMETERS);
		s=getNormalizedSpectrum(s, chargedMz, charge, f, PARAMETERS);
		Charter.launchChart(s);
	}
	
	public void testXCorr() {
		final byte charge=2;
		final double chargedMz=(1329.6335+(charge-1)*MassConstants.protonMass)/charge;
		
		Spectrum s=getSDFHLFGPPGKK();
		
		XCorrCalculator preprocessedSpectrum=new XCorrCalculator(s, chargedMz, charge, PARAMETERS);
		float spectrumFirst=preprocessedSpectrum.score("SDFHLFGPPGKK");
		System.out.println("spectrumFirst xcorr: "+spectrumFirst);
		
		XCorrCalculator preprocessedmodel=new XCorrCalculator("SDFHLFGPPGKK", chargedMz, charge, PARAMETERS);
		float modelFirst=preprocessedmodel.score(s);
		System.out.println("modelFirst xcorr: "+modelFirst);
		
		
		float[] f=XCorrCalculator.normalize(s, s.getPrecursorMZ(), charge, false, PARAMETERS);
		float[] t=XCorrCalculator.getTheoreticalSpectrum("SDFHLFGPPGKK", chargedMz, charge, PARAMETERS);

		System.out.println("-5: "+XCorrCalculator.dotProduct(t, f, -5));
		System.out.println("-4: "+XCorrCalculator.dotProduct(t, f, -4));
		System.out.println("-3: "+XCorrCalculator.dotProduct(t, f, -3));
		System.out.println("-2: "+XCorrCalculator.dotProduct(t, f, -2));
		System.out.println("-1: "+XCorrCalculator.dotProduct(t, f, -1));
		System.out.println(" 0: "+XCorrCalculator.dotProduct(t, f, 0));
		System.out.println(" 1: "+XCorrCalculator.dotProduct(t, f, 1));
		System.out.println(" 2: "+XCorrCalculator.dotProduct(t, f, 2));
		System.out.println(" 3: "+XCorrCalculator.dotProduct(t, f, 3));
		System.out.println(" 4: "+XCorrCalculator.dotProduct(t, f, 4));
		System.out.println(" 5: "+XCorrCalculator.dotProduct(t, f, 5));
		
		float center=XCorrCalculator.dotProduct(t, f, 0);
		float avg=0.0f;
		for (int i=-75; i<75; i++) {
			if (i!=0) {
				avg+=XCorrCalculator.dotProduct(t, f, i);
			}
		}
		avg=avg/150.0f;
		
		float originalCalculation=(center-avg)/1e4f;
		

		float center2=XCorrCalculator.dotProduct(f, t, 0);
		float avg2=0.0f;
		for (int i=-75; i<75; i++) {
			if (i!=0) {
				avg2+=XCorrCalculator.dotProduct(f, t, i);
			}
		}
		avg2=avg2/150.0f;
		
		float originalCalculation2=(center2-avg2)/1e4f;
		assertEquals(originalCalculation, originalCalculation2, 0.05f);
		assertEquals(originalCalculation, modelFirst, 0.05f);
		assertEquals(originalCalculation, spectrumFirst, 0.05f);
		System.out.println(center+"\t"+avg+"\t"+originalCalculation);
		
		//s=getNormalizedSpectrum(s, chargedMz, charge, f, PARAMETERS);
		//Charter.launchChart(s, "model orig:"+originalCalculation2+" spec orig:"+originalCalculation+" model:"+modelFirst+" spec:"+spectrumFirst);
	}

	public static Spectrum getSDFHLFGPPGKK() { 
		// from TRFE_CHICK

		final byte charge=2;
		final double chargedMz=(1329.6335+(charge-1)*MassConstants.protonMass)/charge;
		
		InputStream is=MedianInterpolatorTest.class.getResourceAsStream("/040203_XXX_X1_1_OT_5seq.02.00085.2.dta.txt");
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
	
	static Spectrum getNormalizedSpectrum(final Spectrum s, final double precursorMz, final byte charge, final float[] intensityBins, final SearchParameters params) {
		double massPlusOne=precursorMz*charge-(charge-1)*MassConstants.protonMass;
		// set tolerance to 2x the fragment tolerance of the highest fragment
		float fragmentBinSize=2.0f*(float) params.getFragmentTolerance().getTolerance(massPlusOne);
		double offset;
		if (fragmentBinSize>0.5f) {
			fragmentBinSize=XCorrCalculator.lowResFragmentBinSize; 
			offset=XCorrCalculator.lowResFragmentBinOffset;
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
		for (int i=0; i<intensityBins.length; i++) {
			if (intensityBins[i]!=0.0f) {
				double mass=(i*fragmentBinSize)+offset;
				masses.add(mass);
				intensities.add(intensityBins[i]);
			}
		}
		return getSpectrum(masses.toArray(), intensities.toArray(), tic, scanStartTime, name, mz);
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
