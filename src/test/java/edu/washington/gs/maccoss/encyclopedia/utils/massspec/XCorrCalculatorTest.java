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
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50), DigestionEnzyme.getEnzyme("trypsin"));
	
	public static void main(String[] args) {
		final byte charge=2;
		final double chargedMz=(1329.6335+(charge-1)*MassConstants.protonMass)/charge;
		
		Spectrum s=getSDFHLFGPPGKK();
		float[] f=XCorrCalculator.normalize(s, s.getPrecursorMZ(), charge, false, PARAMETERS);
		f=XCorrCalculator.getTheoreticalSpectrum("SDFHLFGPPGKK", chargedMz, charge, PARAMETERS);
		
		s=getNormalizedSpectrum(s, chargedMz, charge, f, PARAMETERS);
		Charter.launchChart(s);
	}
	
	public void testModel() {
	}

	public static Spectrum getSDFHLFGPPGKK() { 
		// from TRFE_CHICK

		final byte charge=2;
		final double chargedMz=(1329.6335+(charge-1)*MassConstants.protonMass)/charge;
		
		InputStream is=MedianInterpolatorTest.class.getResourceAsStream("/040203_XXX_X1_1_OT_5seq.02.00085.2.dta.txt");
		ArrayList<Peak> peaks=getData(is);
		Pair<double[], float[]> peakArrays=Peak.toArrays(peaks);
		return getSpectrum(peakArrays.x, peakArrays.y, General.sum(peakArrays.y), 0.0f, "SDFHLFGPPGKK", chargedMz); 
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
		
		float fragmentBinSize=2.0f*(float)params.getFragmentTolerance().getTolerance(massPlusOne);
		
		float tic=s.getTIC();
		float scanStartTime=s.getScanStartTime();
		double mz=s.getPrecursorMZ();
		String name=s.getSpectrumName();
		
		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();
		for (int i=0; i<intensityBins.length; i++) {
			if (intensityBins[i]>0.0f) {
				double mass=i*fragmentBinSize;
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
