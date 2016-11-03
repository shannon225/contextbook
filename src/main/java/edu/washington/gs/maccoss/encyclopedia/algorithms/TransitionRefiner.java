package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Optional;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.datastructures.IntRange;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.FloatPair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class TransitionRefiner {
	// minimum threshold to call this peak as worth quantifying
	public static final float quantitativeCorrelationThreshold=0.9f;
	
	// minimum threshold to call this peak as useful for identification purposes
	public static final float identificationCorrelationThreshold=0.75f;
	
	public static void main(String[] args) {
		ArrayList<float[]> chromatograms=new ArrayList<float[]>();
		String[] ionNames=new String[] {"y2", "b3", "b4", "y3", "b5", "y4", "b6", "y5", "y6", "y7", "y8", "y9"};
		float[] y2=new float[] { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 7182.16455078125f, 18434.455078125f, 21684.3671875f, 3613.233642578125f, 8689.09765625f, 12955.7373046875f, 28795.33203125f,
				3359.6435546875f, 7611.09130859375f, 11048.0908203125f, 9528.0302734375f, 12914.23828125f, 8072.17626953125f, 3192.732666015625f, 2322.4375f, 2494.99609375f, 3846.780029296875f,
				3825.619140625f, 2689.070556640625f };
		chromatograms.add(y2);
		float[] b3=new float[] { 23338.361328125f, 16978.677734375f, 26238.6640625f, 28618.11328125f, 47211.97265625f, 60493.10546875f, 85625.6953125f, 154640.59375f, 163637.515625f, 113405.609375f,
				164475.375f, 202257.890625f, 100290.7734375f, 63675.58984375f, 31520.583984375f, 22526.6953125f, 0.0f, 0.0f, 0.0f, 6942.896484375f, 25359.82421875f, 26355.232421875f, 28414.279296875f,
				32256.48046875f, 28046.2421875f };
		chromatograms.add(b3);
		float[] b4=new float[] { 8761.6865234375f, 10261.724609375f, 15003.5693359375f, 15778.66015625f, 13637.501953125f, 9062.185546875f, 8053.20068359375f, 3061.834228515625f, 0.0f,
				429.0523376464844f, 6277.9892578125f, 25745.412109375f, 42210.01171875f, 43571.8984375f, 58320.44140625f, 51980.7578125f, 24411.55078125f, 12893.9833984375f, 4273.0595703125f,
				5722.4169921875f, 3113.877197265625f, 6807.1748046875f, 8823.681640625f, 10008.7529296875f, 20517.623046875f };
		chromatograms.add(b4);
		float[] y3=new float[] { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 15851.9375f, 63487.34375f, 30229.947265625f, 0.0f, 1587.282470703125f, 14519.814453125f, 20986.279296875f, 11382.1982421875f,
				692.5326538085938f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
		chromatograms.add(y3);
		float[] b5=new float[] { 3496.343994140625f, 2400.768798828125f, 3571.334716796875f, 882.3250122070312f, 2112.505615234375f, 0.0f, 0.0f, 0.0f, 8555.8603515625f, 2884.031005859375f,
				6701.763671875f, 7710.79443359375f, 9398.0859375f, 6976.06494140625f, 1426.307373046875f, 0.0f, 1168.6495361328125f, 7716.265625f, 4298.24462890625f, 7817.73779296875f, 0.0f, 0.0f,
				0.0f, 0.0f, 2262.781494140625f };
		chromatograms.add(b5);
		float[] y4=new float[] { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1002.71435546875f, 16844.791015625f, 7548.5625f, 16408.921875f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 0.0f, 0.0f };
		chromatograms.add(y4);
		float[] b6=new float[] { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 4799.2041015625f, 9862.6865234375f, 8904.4296875f,
				2864.54638671875f, 0.0f, 0.0f, 0.0f };
		chromatograms.add(b6);
		float[] y5=new float[] { 0.0f, 0.0f, 18886.91796875f, 63742.88671875f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 22685.7265625f, 46076.71875f, 0.0f, 4507.19189453125f, 0.0f, 0.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
		chromatograms.add(y5);
		float[] y6=new float[] { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 7574.158203125f, 39854.45703125f, 66459.0625f, 26062.421875f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 0.0f, 0.0f };
		chromatograms.add(y6);
		float[] y7=new float[] { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 11383.68359375f, 51658.609375f, 68006.3046875f, 19218.361328125f, 3688.47412109375f, 0.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
		chromatograms.add(y7);
		float[] y8=new float[] { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 14235.25f, 70137.21875f, 104909.390625f, 32442.04296875f, 7820.50927734375f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
		chromatograms.add(y8);
		float[] y9=new float[] { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 14231.3134765625f, 49986.1015625f, 82333.5390625f, 28686.3828125f, 4543.765625f, 0.0f, 0.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
		chromatograms.add(y9);
		
		for (int i=0; i<chromatograms.size(); i++) {
			float[] chromatogram=SkylineSGFilter.paddedSavitzkyGolaySmooth(chromatograms.get(i));
			chromatograms.set(i, chromatogram);
		}

		float[] rts=new float[] { 30.342016220092773f, 30.380882263183594f, 30.42345428466797f, 30.462739944458008f, 30.503339767456055f, 30.543596267700195f, 30.583803176879883f, 30.622554779052734f,
				30.664770126342773f, 30.703386306762695f, 30.74576187133789f, 30.78369140625f, 30.825769424438477f, 30.865869522094727f, 30.906150817871094f, 30.945362091064453f, 30.98526954650879f,
				31.024911880493164f, 31.067047119140625f, 31.105043411254883f, 31.147865295410156f, 31.185548782348633f, 31.228477478027344f, 31.267902374267578f, 31.308513641357422f };
		for (int i=0; i<rts.length; i++) {
			rts[i]=rts[i]*60.0f;
		}

		TDoubleArrayList masses=new TDoubleArrayList();
		int count=0;
		for (int i=0; i<chromatograms.size(); i++) {
			masses.add(count++);
		}
		double[] fragmentMasses=masses.toArray();
		
		TransitionRefinementData data=identifyTransitions("ASVAAQQQEEAR", (byte)2, fragmentMasses, chromatograms, rts, Optional.ofNullable((float[])null), true);
		float[] correlations=data.getCorrelationArray();
		float[] integrations=data.getIntegrationArray();
		for (int i=0; i<integrations.length; i++) {
			System.out.println(ionNames[i]+"\t"+correlations[i]+"\t"+integrations[i]);
		}
		Charter.launchCharts("TITLE", getChartPanels(data));
	}

	public static TransitionRefinementData identifyTransitions(String peptideModSeq, byte precursorCharge, double[] fragmentMasses, ArrayList<float[]> chromatograms, float[] retentionTimes) {
		return identifyTransitions(peptideModSeq, precursorCharge, fragmentMasses, chromatograms, retentionTimes, Optional.ofNullable((float[])null), false);
	}

	public static TransitionRefinementData identifyTransitions(String peptideModSeq, byte precursorCharge, double[] fragmentMasses, ArrayList<float[]> chromatograms, float[] retentionTimes, Optional<float[]> medianChromatogram) {
		return identifyTransitions(peptideModSeq, precursorCharge, fragmentMasses, chromatograms, retentionTimes, medianChromatogram, false);
	}
	static TransitionRefinementData identifyTransitions(String peptideModSeq, byte precursorCharge, double[] fragmentMasses, ArrayList<float[]> chromatograms, float[] retentionTimes, Optional<float[]> maybeMedianChromatogram, boolean plot) {
		if (chromatograms.size()==0) return new TransitionRefinementData(peptideModSeq, precursorCharge, new double[0], chromatograms, new float[0], new float[0], new float[0], new float[0], new Range(retentionTimes[0], retentionTimes[retentionTimes.length-1]));
		
		ArrayList<float[]> normalizedChromatograms;
		int maxIndex;
		IntRange indices;
		float[] medianChromatogram;
		
		if (maybeMedianChromatogram.isPresent()&&maybeMedianChromatogram.get().length>0) {
			// already started with a median chromatogram
			medianChromatogram=maybeMedianChromatogram.get();
			maxIndex=0;
			for (int i=1; i<medianChromatogram.length; i++) {
				if (medianChromatogram[i]>medianChromatogram[maxIndex]) {
					maxIndex=i;
				}
			}
			indices=getIndexRange(medianChromatogram, maxIndex);
			normalizedChromatograms=normalizeAndBackgroundSubtract(chromatograms, indices);
			
		} else {
			// start across the entire width
			normalizedChromatograms=normalize(chromatograms);
			
			// find the maximum point
			medianChromatogram=new float[chromatograms.get(0).length];
			maxIndex=0;
			for (int i=0; i<medianChromatogram.length; i++) {
				TFloatArrayList list=new TFloatArrayList();
				for (float[] chromatogram : normalizedChromatograms) {
					if (chromatogram.length>i) {
						list.add(chromatogram[i]);
					} else {
						list.add(0.0f);
					}
				}
				medianChromatogram[i]=QuickMedian.median(list.toArray());
				if (medianChromatogram[i]>medianChromatogram[maxIndex]) {
					maxIndex=i;
				}
			}
			IntRange initialIndices=getIndexRange(medianChromatogram, maxIndex);

			// then refine on the local area
			normalizedChromatograms=normalizeAndBackgroundSubtract(chromatograms, initialIndices);
			
			// find the maximum point
			medianChromatogram=new float[chromatograms.get(0).length];
			maxIndex=0;
			for (int i=0; i<medianChromatogram.length; i++) {
				TFloatArrayList list=new TFloatArrayList();
				for (float[] chromatogram : normalizedChromatograms) {
					if (chromatogram.length>i) {
						list.add(chromatogram[i]);
					} else {
						list.add(0.0f);
					}
				}
				medianChromatogram[i]=QuickMedian.median(list.toArray());
				if (medianChromatogram[i]>medianChromatogram[maxIndex]) {
					maxIndex=i;
				}
			}
			indices=getIndexRange(medianChromatogram, maxIndex);
		}

		float medianMean=General.mean(medianChromatogram, indices.getStart(), indices.getStop());
		float[] correlationArray=new float[normalizedChromatograms.size()];
		float[] integrationArray=new float[correlationArray.length];
		float[] backgroundArray=new float[correlationArray.length];
		for (int i=0; i<normalizedChromatograms.size(); i++) {
			float[] normalizedChromatogram=normalizedChromatograms.get(i);
			float correlation=calculateCorrelation(medianMean, indices, medianChromatogram, normalizedChromatogram);
			correlationArray[i]=correlation;

			FloatPair intensity=integrate(indices, retentionTimes, chromatograms.get(i));

			integrationArray[i]=intensity.getOne();
			backgroundArray[i]=intensity.getTwo();
			
			// calculate trapezoidal background area
			integrationArray[i]=integrationArray[i]-backgroundArray[i];
		}

		Range range=new Range(retentionTimes[indices.getStart()], retentionTimes[indices.getStop()]);
		if (plot) {
			HashMap<String, ChartPanel> panels=new HashMap<String, ChartPanel>();
			panels.put("unnormalized", getChart(chromatograms, correlationArray, retentionTimes, range));
			panels.put("unnormalized_uncolored", getChart(chromatograms, new float[correlationArray.length], retentionTimes, range));
			panels.put("normalized", getChart(normalizedChromatograms, correlationArray, retentionTimes, range));
			panels.put("median", Charter.getChart("scan", "intensity", false, toXYTrace(medianChromatogram, retentionTimes, "median", null, null)));
			
			ArrayList<XYTrace> traces=getTraces(normalizedChromatograms, correlationArray, retentionTimes, range);
			traces.add(0, toXYTrace(medianChromatogram, retentionTimes, "median", Color.black, null, GraphType.dashedline, 6.0f));
			panels.put("traces", Charter.getChart("scan", "intensity", false, traces.toArray(new XYTrace[traces.size()])));
			
			Charter.launchCharts(peptideModSeq+" chart", panels);
		}
		
		return new TransitionRefinementData(peptideModSeq, precursorCharge, fragmentMasses, chromatograms, correlationArray, integrationArray, backgroundArray, medianChromatogram, range);
	}
	
	public static int[] numberOfCoelutingIons(double[] targetMasses, double[] allIons, ArrayList<Spectrum> stripes, int halfPeakWidthInScans, MassTolerance fragmentTolerance) {
		float[][] targetIntensityArray=new float[stripes.size()][];
		float[][] allIntensityArray=new float[stripes.size()][];
		for (int rtIndex=0; rtIndex<stripes.size(); rtIndex++) {
			Spectrum spectrum=stripes.get(rtIndex);
			float[] allIntegratedIntensities=fragmentTolerance.getIntegratedIntensities(spectrum.getMassArray(), spectrum.getIntensityArray(), allIons);
			allIntensityArray[rtIndex]=allIntegratedIntensities;
			float[] targetIntegratedIntensities=fragmentTolerance.getIntegratedIntensities(spectrum.getMassArray(), spectrum.getIntensityArray(), targetMasses);
			targetIntensityArray[rtIndex]=targetIntegratedIntensities;
		}

		float[][] targetChromatograms=extractChromatograms(targetIntensityArray);
		float[][] allChromatograms=extractChromatograms(allIntensityArray);
		
		int[] complementaryIons=new int[stripes.size()];
		for (int rtIndex=halfPeakWidthInScans; rtIndex<targetIntensityArray.length-halfPeakWidthInScans; rtIndex++) {
			IntRange indexRange=new IntRange(rtIndex-halfPeakWidthInScans, rtIndex+halfPeakWidthInScans);
			
			float[][] targetNormalized=new float[targetChromatograms.length][];
			for (int ionIndex=0; ionIndex<targetNormalized.length; ionIndex++) {
				targetNormalized[ionIndex]=General.normalize(General.extract(targetChromatograms[ionIndex], indexRange));
			}
			
			// calculate median
			float[] median=new float[indexRange.getRange()];
			for (int localRtIndex=0; localRtIndex<median.length; localRtIndex++) {
				TFloatArrayList ions=new TFloatArrayList();
				for (int ionIndex=0; ionIndex<targetNormalized.length; ionIndex++) {
					if (targetNormalized[ionIndex][localRtIndex]>0.0f) {
						ions.add(targetNormalized[ionIndex][localRtIndex]);
					}
				}
				median[localRtIndex]=ions.size()==0?0.0f:QuickMedian.median(ions.toArray());
				//medianMap.put(stripes.get(rtIndex+localRtIndex-halfPeakWidthInScans).getScanStartTime(), median[localRtIndex]);
			}
			//traces.add(new XYTrace(medianMap, GraphType.line, Float.toString(stripes.get(rtIndex).getScanStartTime())));
			
			// normalize and score all chromatograms
			IntRange completeLocalRange=new IntRange(0, median.length-1);
			float medianMean=General.mean(median);
			for (int ionIndex=0; ionIndex<targetNormalized.length; ionIndex++) {
				float[] normalizedChromatogram=General.normalize(General.extract(allChromatograms[ionIndex], indexRange));
				float correlation=TransitionRefiner.calculateCorrelation(medianMean, completeLocalRange, median, normalizedChromatogram);
				if (correlation>=TransitionRefiner.identificationCorrelationThreshold) {
					complementaryIons[rtIndex]++;
				}
			}
		}

		int[] maxIons=new int[complementaryIons.length];
		for (int rtIndex=halfPeakWidthInScans; rtIndex<targetIntensityArray.length-halfPeakWidthInScans; rtIndex++) {
			IntRange indexRange=new IntRange(rtIndex-halfPeakWidthInScans, rtIndex+halfPeakWidthInScans);
			for (int localRtIndex=0; localRtIndex<indexRange.getRange(); localRtIndex++) {
				maxIons[rtIndex+localRtIndex-halfPeakWidthInScans]=General.max(General.extract(complementaryIons, indexRange));
			}
		}
		
		return maxIons;
	}
	
	static float[][] extractChromatograms(float[][] intensities) {
		float[][] chromatograms=new float[intensities[0].length][];
		for (int ionIndex=0; ionIndex<chromatograms.length; ionIndex++) {
			float[] chromatogram=new float[intensities.length];
			for (int rtIndex=0; rtIndex<chromatogram.length; rtIndex++) {
				chromatogram[rtIndex]=intensities[rtIndex][ionIndex];
			}
			chromatograms[ionIndex]=chromatogram;
		}
		return chromatograms;
	}

	public static float calculateCorrelation(float medianMean, IntRange indices, float[] medianChromatogram, float[] normalizedChromatogram) {
		float fragmentMean=General.mean(normalizedChromatogram, indices.getStart(), indices.getStop());
		
		float medianDeltaSquareSum=0.0f;
		float fragmentDeltaSquareSum=0.0f;
		float deltaProductSum=0.0f;
		for (int j=indices.getStart(); j<=indices.getStop(); j++) {
			float deltaMedian=medianChromatogram[j]-medianMean;
			float deltaFragment=normalizedChromatogram[j]-fragmentMean;
			medianDeltaSquareSum+=deltaMedian*deltaMedian;
			fragmentDeltaSquareSum+=deltaFragment*deltaFragment;
			deltaProductSum+=deltaMedian*deltaFragment;
		}
		float correlation;
		// calculate correlation
		if (fragmentDeltaSquareSum==0.0f) {
			correlation=Float.MIN_VALUE;
		} else if (medianDeltaSquareSum==0.0f) {
			correlation=Float.MIN_VALUE;
		} else {
			float denominator=(float)Math.sqrt(medianDeltaSquareSum*fragmentDeltaSquareSum);
			correlation=deltaProductSum/denominator;
			if (correlation>1.0f) {
				correlation=1.0f; // there can be minor floating point errors in the sqrt
			}
		}
		return correlation;
	}

	/**
	 * 
	 * @param indices
	 * @param retentionTimes
	 * @param chromatogram
	 * @return FloatPair<Intensity,Background>
	 */
	public static FloatPair integrate(IntRange indices, float[] retentionTimes, float[] chromatogram) {
		float background=0.0f;
		float integration=0.0f;
		float baseBackground=Math.min(chromatogram[indices.getStart()], chromatogram[indices.getStop()]);
		
		for (int j=indices.getStart()+1; j<=indices.getStop(); j++) {
			// calculate trapezoid area like Reimann sums
			float trapezoid=(retentionTimes[j]-retentionTimes[j-1])*(chromatogram[j-1]+chromatogram[j])/2.0f;
			integration+=trapezoid;

			// Skyline computes rectangular background based on minimum of the boundaries. If the signal dips below the background, then calculate the trapezoid
			float lowerBackground=Math.min(baseBackground, chromatogram[j-1]);
			float upperBackground=Math.min(baseBackground, chromatogram[j]);
			float trapezoidBackground=(retentionTimes[j]-retentionTimes[j-1])*(lowerBackground+upperBackground)/2.0f;
			background+=trapezoidBackground;
		}
		FloatPair intensity=new FloatPair(integration, background);
		return intensity;
	}

	private static IntRange getIndexRange(float[] medianChromatogram, int maxIndex) {
		float threshold=medianChromatogram[maxIndex]*0.01f; // 1% of max
		
		// left of center (decreasing index)
		int increasing=0;
		int firstData=maxIndex;
		for (int i=maxIndex-1; i>=0; i--) {
			// navigate down the slope of the peak:
			// count the number of consecutive uphill points (count the aggregate, so +1 for increasing, -1 for decreasing)
			if (medianChromatogram[i]>medianChromatogram[firstData]) {
				increasing++;
			} else if (increasing>0) {
				increasing--;
			}
			
			// create peak boundary if we've seen 3 or more consecutively increasing points and we're less than 50% of the max
			if (increasing>2&&medianChromatogram[maxIndex]/2.0f>medianChromatogram[firstData]) {
				break;
			}
			
			// if we're lower than the previous local minimum, set the new minimum
			if (medianChromatogram[i]<medianChromatogram[firstData]) {
				firstData=i;
			}
			
			// create peak boundary if the local minimum is less than 1%
			if (medianChromatogram[firstData]<threshold) {
				break;
			}
		}
		
		// right of center (increasing index)
		increasing=0;
		int lastData=maxIndex;
		for (int i=maxIndex+1; i<medianChromatogram.length; i++) {
			// navigate down the slope of the peak:
			// count the number of consecutive uphill points (count the aggregate, so +1 for increasing, -1 for decreasing)
			if (medianChromatogram[i]>medianChromatogram[lastData]) {
				increasing++;
			} else if (increasing>0) {
				increasing--;
			}
			
			// create peak boundary if we've seen 3 or more consecutively increasing points and we're less than 50% of the max
			if (increasing>2&&medianChromatogram[maxIndex]/2.0f>medianChromatogram[lastData]) {
				break;
			}

			// if we're lower than the previous local minimum, set the new minimum
			if (medianChromatogram[i]<medianChromatogram[lastData]) {
				lastData=i;
			}
			
			// create peak boundary if the local minimum is less than 1%
			if (medianChromatogram[lastData]<threshold) {
				break;
			}
		}
		
		int startIndex=firstData<=0?0:firstData;
		int stopIndex=lastData>=medianChromatogram.length-1?medianChromatogram.length-1:lastData;
		IntRange indices=new IntRange(startIndex, stopIndex);
		return indices;
	}
	
	public static HashMap<String, ChartPanel> getChartPanels(TransitionRefinementData data) {
		HashMap<String, ChartPanel> panels=new HashMap<String, ChartPanel>();
		float[] rts=null;
		if (data!=null&&data.getRtArray()!=null&&data.getRtArray().isPresent()) {
			rts=data.getRtArray().get();
		}
		panels.put("unnormalized", getChart(data.getChromatograms(), data.getCorrelationArray(), rts, data.getRange()));
		panels.put("median", Charter.getChart("Retention Time (min)", "intensity", false, toXYTrace(data.getMedianChromatogram(), rts, "median", null, data.getRange())));
		return panels;
	}

	private static ChartPanel getChart(ArrayList<float[]> chromatograms, float[] correlationArray, float[] rts, Range rtRange) {
		ArrayList<XYTrace> xytraces=getTraces(chromatograms, correlationArray, rts, rtRange);
		
		return tracesToChart(xytraces);
	}

	public static ChartPanel tracesToChart(ArrayList<XYTrace> xytraces) {
		ChartPanel panel=Charter.getChart("Retention Time (min)", "intensity", false, xytraces.toArray(new XYTrace[xytraces.size()]));
		return panel;
	}

	private static ArrayList<XYTrace> getTraces(ArrayList<float[]> chromatograms, float[] correlationArray, float[] rts, Range rtRange) {
		ArrayList<XYTrace> xytraces=new ArrayList<XYTrace>();
		for (int i=0; i<chromatograms.size(); i++) {
			float[] fs=chromatograms.get(i);
			
			Color c;
			if (correlationArray[i]>quantitativeCorrelationThreshold) {
				c=new Color(0, 205, 0);
			} else if (correlationArray[i]>identificationCorrelationThreshold) {
				c=new Color(255, 215, 0);
			} else if (correlationArray[i]==0.0f) {
				c=Color.gray;
			} else {
				c=Color.red;
			}

			xytraces.add(toXYTrace(fs, rts, ""+i, c, rtRange));

			if (rtRange!=null) {
				xytraces.add(toXYTrace(fs, rts, ""+i, Color.gray, null));
			}
		}
		return xytraces;
	}
	public static XYTraceInterface toBoundaries(float f, String name) {
		ArrayList<XYPoint> points=new ArrayList<XYPoint>();
		points.add(new XYPoint(f, -0.1f));
		points.add(new XYPoint(f, 0.1f));
		return new XYTrace(points, GraphType.line, name);
	}

	public static XYTrace toXYTrace(float[] fs, float[] rts, String name, Color color, Range rtRange) {
		GraphType graphtype=GraphType.line;
		float thickness=3.0f;
		return toXYTrace(fs, rts, name, color, rtRange, graphtype, thickness);
	}

	public static XYTrace toXYTrace(float[] fs, float[] rts, String name, Color color, Range rtRange, GraphType graphtype, float thickness) {
		ArrayList<XYPoint> points=new ArrayList<XYPoint>();
		for (int j=0; j<fs.length; j++) {
			if (rts==null) {
				XYPoint point=new XYPoint(j, fs[j]);
				points.add(point);
			} else {
				if (rtRange!=null&&!rtRange.contains(rts[j])) continue;
				
				XYPoint point=new XYPoint(rts[j]/60f, fs[j]);
				points.add(point);
			}
		}
		XYTrace trace=new XYTrace(points, graphtype, name, color, thickness);
		return trace;
	}

	public static ArrayList<float[]> normalize(ArrayList<float[]> chromatograms) {
		ArrayList<float[]> normalizedChromatograms=new ArrayList<float[]>();
		for (float[] fs : chromatograms) {
			normalizedChromatograms.add(General.normalize(fs));
		}
		return normalizedChromatograms;
	}

	public static ArrayList<float[]> normalizeAndBackgroundSubtract(ArrayList<float[]> chromatograms, IntRange range) {
		ArrayList<float[]> normalizedChromatograms=new ArrayList<float[]>();
		for (float[] fs : chromatograms) {
			// TODO CONSIDER PUTTING BACKGROUND SUBTRACTION INTO CORRELATION! (but not this way)
			//normalizedChromatograms.add(General.normalizeAndBackgroundSubtract(fs, range));
			normalizedChromatograms.add(General.normalize(fs, range));
		}
		return normalizedChromatograms;
	}
}
