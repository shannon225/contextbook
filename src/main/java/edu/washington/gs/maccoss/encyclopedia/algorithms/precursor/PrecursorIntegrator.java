package edu.washington.gs.maccoss.encyclopedia.algorithms.precursor;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntRange;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeakTrace;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.RetentionTimeBoundary;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Ion;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PrecursorIon;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Correlation;
import edu.washington.gs.maccoss.encyclopedia.utils.math.FloatPair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TByteObjectHashMap;
import gnu.trove.procedure.TByteObjectProcedure;

public class PrecursorIntegrator {
	private static final float REQUIRED_ION_PERCENTAGE_OF_DISTRIBUTION = 0.5f;
	
	@SuppressWarnings("unchecked")
	public static ArrayList<TransitionRefinementData> integratePeptide(String peptideModSeq, float[] expectedIsotopicDistribution, float[] psmRTsInSec, Optional<Float> ionMobility, ArrayList<PeakTrace<PrecursorIon>> traces, SearchParameters params) {
		TByteObjectHashMap<PeakTrace<PrecursorIon>[]> tracesByCharge = getTracesByCharge(traces);

		ArrayList<PeakTrace<PrecursorIon>> distributionCorrectedTraces = new ArrayList<>();
		// for each charge state
		tracesByCharge.forEachEntry(new TByteObjectProcedure<PeakTrace<PrecursorIon>[]>() {
			@Override
			public boolean execute(byte charge, PeakTrace<PrecursorIon>[] traces) {
				// remove the -1 trace for quant (it is only used to confirm detection since there shouldn't be any -1 signal)
				PeakTrace<PrecursorIon>[] tracesCopy=new PeakTrace[traces.length-1];
				System.arraycopy(traces, 1, tracesCopy, 0, tracesCopy.length);
				
				int numPoints=tracesCopy[0].size();
				
				TFloatArrayList[] pureTraces=new TFloatArrayList[tracesCopy.length];
				for (int i = 0; i < pureTraces.length; i++) {
					pureTraces[i]=new TFloatArrayList();
				} 
				
				for (int t = 0; t < numPoints; t++) { // for every RT time point (t)
					float[] intensities=new float[tracesCopy.length];
					for (int i=0; i<tracesCopy.length; i++) {
						intensities[i]=tracesCopy[i].getIntensity()[t];
					}
					
					float minScaler=Float.MAX_VALUE;
					for (int i=0; i<intensities.length; i++) {
						if (expectedIsotopicDistribution[i]>REQUIRED_ION_PERCENTAGE_OF_DISTRIBUTION) {
							float scaler=intensities[i]/expectedIsotopicDistribution[i];
							if (minScaler>scaler) {
								minScaler=scaler;
							}
							//System.out.println(charge+","+traces[i].getRt()[t]+" = "+intensities[i]+" vs "+(expectedIsotopicDistribution[i]*minScaler)+", ("+expectedIsotopicDistribution[i]+", "+scaler+")");
						}
					}

					for (int i=0; i<tracesCopy.length; i++) {
						// can't ever have more intensity than recorded
						pureTraces[i].add(Math.min(intensities[i], expectedIsotopicDistribution[i]*minScaler));
					}
				}

				for (int i = 0; i < pureTraces.length; i++) {
					distributionCorrectedTraces.add(tracesCopy[i].updateIntensity(pureTraces[i].toArray()));
				}
				return true;
			}
		});

		Pair<ArrayList<PeakTrace<PrecursorIon>>, ArrayList<PeakTrace<PrecursorIon>>> smoothedAndNormalized = smoothAndNormalize(distributionCorrectedTraces, true, psmRTsInSec, params);

		@SuppressWarnings("rawtypes")
		ArrayList<PeakTrace> genericSmoothedAndNormalized=generify(smoothedAndNormalized.y);
		PeakTrace<String> median=PeakTrace.getMedianTrace(genericSmoothedAndNormalized);
		RetentionTimeBoundary boundary = getBoundaries(median, psmRTsInSec, params);
		
		if (0.0f==General.max(median.getIntensity())) {
			ArrayList<PeakTrace<PrecursorIon>> originalTraces=new ArrayList<>();
			for (PeakTrace<PrecursorIon> peakTrace : traces) {
				if (peakTrace.getIon().getIsotope()>=0&&peakTrace.getIon().getIsotope()<=1) { // work from monoisotopic and +1 only
					originalTraces.add(peakTrace);
				}
			} 
			smoothedAndNormalized = smoothAndNormalize(originalTraces, false, psmRTsInSec, params);
			genericSmoothedAndNormalized=generify(smoothedAndNormalized.y);
			median=PeakTrace.getMedianTrace(genericSmoothedAndNormalized); // perhaps sub with something that requires starting from a certain RT?
			boundary = getBoundaries(median, psmRTsInSec, params);
		}
		
		final PeakTrace<String> finalMedian=PeakTrace.getMedianTrace(genericSmoothedAndNormalized, boundary);
		final float[] rtArray=finalMedian.getRt();
		final float[] medianIntensityArray=finalMedian.getIntensity();
		final RetentionTimeBoundary finalBoundary=boundary;
		
		TByteObjectHashMap<PeakTrace<PrecursorIon>[]> smoothedTracesByCharge = getTracesByCharge(smoothedAndNormalized.x);

		final ArrayList<TransitionRefinementData> data=new ArrayList<>();
		smoothedTracesByCharge.forEachEntry(new TByteObjectProcedure<PeakTrace<PrecursorIon>[]>() {
			@Override
			public boolean execute(byte charge, PeakTrace<PrecursorIon>[] smoothedTraces) {
				ArrayList<Ion> ions=new ArrayList<>(); // FRAGMENTION
				TDoubleArrayList fragmentMasses=new TDoubleArrayList();
				TFloatArrayList deltaMassArray=new TFloatArrayList();
				TFloatArrayList correlationArray=new TFloatArrayList();
				TFloatArrayList integrationArray=new TFloatArrayList();
				TFloatArrayList backgroundArray=new TFloatArrayList();
				ArrayList<float[]> chromatograms=new ArrayList<>();
				for (int i = 0; i < smoothedTraces.length; i++) {
					if (smoothedTraces[i]!=null) {
						ions.add(smoothedTraces[i].getIon());
						deltaMassArray.add(0.0f); // FIXME add delta masses!
						fragmentMasses.add(smoothedTraces[i].getIon().getMass());
						float correlation=smoothedTraces[i].getCorrelation(finalMedian, finalBoundary);
						correlationArray.add(correlation);
						FloatPair intensity=smoothedTraces[i].integrateIntensityAndBackground(finalBoundary);
						integrationArray.add(intensity.getOne());
						backgroundArray.add(intensity.getTwo());
						chromatograms.add(smoothedTraces[i].getIntensity());
					}
				}
				boolean[] quantitativeIonsArray=new boolean[correlationArray.size()];
				for (int i = 0; i < quantitativeIonsArray.length; i++) {
					quantitativeIonsArray[i]=correlationArray.get(i)>TransitionRefiner.quantitativeCorrelationThreshold;
				}
				TransitionRefinementData trd=new TransitionRefinementData(peptideModSeq, charge, ions.toArray(new Ion[ions.size()]), chromatograms, 
						correlationArray.toArray(), quantitativeIonsArray, integrationArray.toArray(), backgroundArray.toArray(), medianIntensityArray, ionMobility, finalBoundary, 
						deltaMassArray.toArray(), fragmentMasses.toArray(), integrationArray.toArray(), rtArray, null, null, 0.0f, 
						params.getAAConstants());
				data.add(trd);
				
				return true;
			}
		});
		return data;
		
//		final ArrayList<IntegrationScores> allScores=new ArrayList<>();
//		smoothedTracesByCharge.forEachEntry(new TByteObjectProcedure<PeakTrace<PrecursorIon>[]>() {
//			@Override
//			public boolean execute(byte charge, PeakTrace<PrecursorIon>[] smoothedTraces) {
//				TFloatArrayList correlationList=new TFloatArrayList();
//				TFloatArrayList apexRTList=new TFloatArrayList();
//				TFloatArrayList intensityList=new TFloatArrayList();
//				for (PeakTrace<PrecursorIon> trace : smoothedTraces) {
//					if (trace==null) continue; // for -1 trace
//					
//					float correlation=trace.getCorrelation(finalMedian, finalBoundary);
//					correlationList.add(correlation);
//					float apexRT = trace.getApexRT(finalBoundary);
//					apexRTList.add(apexRT);
//					float intensity=trace.integrateIntensity(finalBoundary, true);
//					intensityList.add(intensity);
//					//System.out.println(trace.getIon()+") Corr:"+correlation+", RT:"+apexRT+", Int:"+intensity+" --> "+trace.multiply(2).getCorrelation(finalMedian, boundary));
//				}
//				float averagePeakCorrelation=General.mean(correlationList.toArray());
//				float meanApexRT=General.mean(apexRTList.toArray());
//				float totalIntensity=General.sum(intensityList.toArray());
//				
//				PeakTrace<PrecursorIon>[] originalTraces=tracesByCharge.get(charge);
//				float[] originalIntensities=new float[originalTraces.length];
//				float[] expectedDistribution=new float[originalIntensities.length];
//				for (int i = 0; i < originalTraces.length; i++) {
//					byte isotope=originalTraces[i].getIon().getIsotope();
//					originalIntensities[i]=originalTraces[i].integrateIntensity(finalBoundary, true);
//					if (isotope>=0) {
//						expectedDistribution[i]=expectedIsotopicDistribution[isotope];
//					}
//				}
//				float isotopicRatioCorrelation=Correlation.getPearsons(originalIntensities, expectedDistribution);
//
//				IntegrationScores scores=new IntegrationScores(smoothedTraces[1].getIon(), totalIntensity, meanApexRT, averagePeakCorrelation, isotopicRatioCorrelation, finalBoundary, finalMedian);
//				allScores.add(scores);
//				
//				return true;
//			}
//		});
//		
//		return allScores;
	}

	private static Pair<ArrayList<PeakTrace<PrecursorIon>>, ArrayList<PeakTrace<PrecursorIon>>> smoothAndNormalize(
			ArrayList<PeakTrace<PrecursorIon>> distributionCorrectedTraces, boolean useSmoothed, float[] psmRTsInSec, SearchParameters params) {
		ArrayList<PeakTrace<PrecursorIon>> smoothedList=new ArrayList<>();
		ArrayList<PeakTrace<PrecursorIon>> smoothedAndNormalizedList=new ArrayList<>();
		for (PeakTrace<PrecursorIon> trace : distributionCorrectedTraces) {
			if (((PrecursorIon)trace.getIon()).getIsotope()>=0) {
				PeakTrace<PrecursorIon> smooth;
				if (useSmoothed) {
					smooth=trace.smooth();
				} else {
					smooth=trace;
				}
				RetentionTimeBoundary b=getBoundaries(smooth, psmRTsInSec, params);
				float apexIntensity=smooth.getIntensity(smooth.getApexRT(b));
				smoothedList.add(smooth);
				smoothedAndNormalizedList.add(smooth.multiply(1f/apexIntensity));
			}
		}
		
		Pair<ArrayList<PeakTrace<PrecursorIon>>,ArrayList<PeakTrace<PrecursorIon>>> smoothedAndNormalized=new Pair<ArrayList<PeakTrace<PrecursorIon>>, ArrayList<PeakTrace<PrecursorIon>>>(smoothedList, smoothedAndNormalizedList);
		return smoothedAndNormalized;
	}

	@SuppressWarnings("unchecked")
	private static TByteObjectHashMap<PeakTrace<PrecursorIon>[]> getTracesByCharge(ArrayList<PeakTrace<PrecursorIon>> traces) {
		TByteObjectHashMap<PeakTrace<PrecursorIon>[]> tracesByCharge=new TByteObjectHashMap<>();
		for (PeakTrace<PrecursorIon> trace : traces) {
			PrecursorIon ion=trace.getIon();
			PeakTrace<PrecursorIon>[] array=tracesByCharge.get(ion.getCharge());
			if (array==null) {
				array=new PeakTrace[IntegratedPeptide.integratedIsotopes.length];
				tracesByCharge.put(ion.getCharge(), array);
			}
			array[ion.getIsotope()+1]=trace; // keep the -1 isotope
		}
		return tracesByCharge;
	}
	
	@SuppressWarnings("rawtypes")
	private static ArrayList<PeakTrace> generify(ArrayList<PeakTrace<PrecursorIon>> traces) {
		return new ArrayList<PeakTrace>(traces);
	}

	private static RetentionTimeBoundary getBoundaries(@SuppressWarnings("rawtypes") PeakTrace peak, float[] psmRTsInSec, SearchParameters params) {
		float[] rts=peak.getRt();
		float[] intensities=peak.getIntensity();
		
		IntRange range=TransitionRefiner.getIndexRange(rts, intensities, params.getExpectedPeakWidth(), params.getMinNumIntegratedRTPoints());
		
		// require that the boundary contains the PSMs
		RetentionTimeBoundary boundaryFromPSMs = RetentionTimeBoundary.getMedianBoundaries(psmRTsInSec);
		float start = Math.min(boundaryFromPSMs.getMinRT(), rts[range.getStart()]);
		float stop = Math.max(boundaryFromPSMs.getMaxRT(), rts[range.getStop()]);
		Range boundaryRange=new Range(start, stop);
		
		int maxIndex=0;
		for (int i = 1; i < intensities.length; i++) {
			if (boundaryRange.contains(rts[i])) {
				if (intensities[i]>intensities[maxIndex]) {
					maxIndex=i;
				}
			}
		}
		
		RetentionTimeBoundary boundary=new RetentionTimeBoundary(start, rts[maxIndex], stop);
		return boundary;
	}
}