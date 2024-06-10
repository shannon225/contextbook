package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Optional;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntRange;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.FloatPair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class TransitionRefiner {
	private static final float NO_AVAILABLE_INTENSITY_MARKER = -1.0f;
	public static final float PERCENT_RT_UPHILL_PEAKWIDTH = 0.2f;
	public static final float MINIMUM_THRESHOLD_PERCENTAGE = 0.01f;

	// minimum threshold to call this peak as worth quantifying
	public static final float quantitativeCorrelationThreshold=0.9f;
	
	// minimum threshold to call this peak as useful for identification purposes
	public static final float identificationCorrelationThreshold=0.75f;
	
	public static void main(String[] args) {
		ArrayList<float[]> chromatograms=new ArrayList<float[]>();
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

		String[] ionNames=new String[] {"y2", "b3", "b4", "y3", "b5", "y4", "b6", "y5", "y6", "y7", "y8", "y9"};
		ArrayList<FragmentIon> ions=new ArrayList<FragmentIon>();
		ions.add(new FragmentIon(2, (byte)2, IonType.y));
		ions.add(new FragmentIon(3, (byte)3, IonType.b));
		ions.add(new FragmentIon(4, (byte)4, IonType.b));
		ions.add(new FragmentIon(3, (byte)3, IonType.y));
		ions.add(new FragmentIon(5, (byte)5, IonType.b));
		ions.add(new FragmentIon(4, (byte)4, IonType.y));
		ions.add(new FragmentIon(6, (byte)6, IonType.b));
		ions.add(new FragmentIon(5, (byte)5, IonType.y));
		ions.add(new FragmentIon(6, (byte)6, IonType.y));
		ions.add(new FragmentIon(7, (byte)7, IonType.y));
		ions.add(new FragmentIon(8, (byte)8, IonType.y));
		ions.add(new FragmentIon(9, (byte)9, IonType.y));
		FragmentIon[] fragmentMasses=ions.toArray(new FragmentIon[ions.size()]);
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		
		TransitionRefinementData data=identifyTransitions("ASVAAQQQEEAR", (byte)2, 30.74576187133789f, fragmentMasses, chromatograms, rts, Optional.ofNullable((float[])null), true, params);
		float[] correlations=data.getCorrelationArray();
		float[] integrations=data.getIntegrationArray();
		for (int i=0; i<integrations.length; i++) {
			System.out.println(ionNames[i]+"\t"+correlations[i]+"\t"+integrations[i]);
		}
		Charter.launchCharts("TITLE", getChartPanels(data));
	}
	
	public static TransitionRefinementData identifyTransitionsFromRTRange(String peptideModSeq, byte precursorCharge, float retentionTimeInSec, FragmentIon[] fragmentMasses, ArrayList<float[]> chromatograms, float[] retentionTimes, Range rtRange, SearchParameters params) {
		if (chromatograms.size()==0) return new TransitionRefinementData(peptideModSeq, precursorCharge, new FragmentIon[0], chromatograms, new float[0], new boolean[0], new float[0], new float[0], new float[0], new Range(retentionTimes[0], retentionTimes[retentionTimes.length-1]), params.getAAConstants());
		
		MedianChromatogramData medianData = extractMedianChromatogramFromFixedBoundaries(retentionTimeInSec, chromatograms, retentionTimes, rtRange);

		float medianMean=General.mean(medianData.getMedianChromatogram(), medianData.getIndices().getStart(), medianData.getIndices().getStop());
		float[] correlationArray=new float[medianData.getNormalizedChromatograms().size()];
		boolean[] quantitativeIonsArray=new boolean[correlationArray.length];
		
		float[] integrationArray=new float[correlationArray.length];
		float[] backgroundArray=new float[correlationArray.length];
		for (int i=0; i<medianData.getNormalizedChromatograms().size(); i++) {
			float[] normalizedChromatogram=medianData.getNormalizedChromatograms().get(i);
			float correlation=calculateCorrelation(medianMean, medianData.getIndices(), medianData.getMedianChromatogram(), normalizedChromatogram);
			correlationArray[i]=correlation;
			quantitativeIonsArray[i]=correlation>=quantitativeCorrelationThreshold;

			FloatPair intensity=integrate(medianData.getIndices(), retentionTimes, chromatograms.get(i));

			integrationArray[i]=intensity.getOne();
			backgroundArray[i]=intensity.getTwo();
			
			// calculate trapezoidal background area
			integrationArray[i]=integrationArray[i]-backgroundArray[i];
		}

		Range range=new Range(retentionTimes[medianData.getIndices().getStart()], retentionTimes[medianData.getIndices().getStop()]);
		
		return new TransitionRefinementData(peptideModSeq, precursorCharge, fragmentMasses, chromatograms, correlationArray, quantitativeIonsArray, integrationArray, backgroundArray, medianData.getMedianChromatogram(), range, params.getAAConstants());
	}

	public static TransitionRefinementData identifyTransitions(String peptideModSeq, byte precursorCharge, float retentionTimeInSec, FragmentIon[] fragmentMasses, ArrayList<float[]> chromatograms, float[] retentionTimes, boolean wasInferred, SearchParameters params) {
		return identifyTransitions(peptideModSeq, precursorCharge, retentionTimeInSec, fragmentMasses, chromatograms, retentionTimes, Optional.ofNullable((float[])null), wasInferred, false, params);
	}

	public static TransitionRefinementData identifyTransitions(String peptideModSeq, byte precursorCharge, float retentionTimeInSec, FragmentIon[] fragmentMasses, ArrayList<float[]> chromatograms, float[] retentionTimes, Optional<float[]> medianChromatogram, boolean wasInferred, SearchParameters params) {
		return identifyTransitions(peptideModSeq, precursorCharge, retentionTimeInSec, fragmentMasses, chromatograms, retentionTimes, medianChromatogram, wasInferred, false, params);
	}

	public static TransitionRefinementData identifyTransitions(String peptideModSeq, byte precursorCharge, float retentionTimeInSec, FragmentIon[] fragmentMasses, ArrayList<float[]> chromatograms, float[] retentionTimes, Optional<float[]> maybeMedianChromatogram, boolean wasInferred, boolean plot, SearchParameters params) {
		if (chromatograms.size()==0) return new TransitionRefinementData(peptideModSeq, precursorCharge, new FragmentIon[0], chromatograms, new float[0], new boolean[0], new float[0], new float[0], new float[0], new Range(retentionTimes[0], retentionTimes[retentionTimes.length-1]), params.getAAConstants());
		
		boolean adjustPeakBoundaries=wasInferred&&params.adjustInferredRTBoundaries();
		MedianChromatogramData medianData = extractMedianChromatogram(retentionTimeInSec, chromatograms, retentionTimes, maybeMedianChromatogram, adjustPeakBoundaries, params.getMinNumIntegratedRTPoints(), params.getExpectedPeakWidth());

		float medianMean=General.mean(medianData.getMedianChromatogram(), medianData.getIndices().getStart(), medianData.getIndices().getStop());
		float[] correlationArray=new float[medianData.getNormalizedChromatograms().size()];
		boolean[] quantitativeIonsArray=new boolean[correlationArray.length];
		
		float[] integrationArray=new float[correlationArray.length];
		float[] backgroundArray=new float[correlationArray.length];
		for (int i=0; i<medianData.getNormalizedChromatograms().size(); i++) {
			float[] normalizedChromatogram=medianData.getNormalizedChromatograms().get(i);
			float correlation=calculateCorrelation(medianMean, medianData.getIndices(), medianData.getMedianChromatogram(), normalizedChromatogram);
			correlationArray[i]=correlation;
			quantitativeIonsArray[i]=correlation>=quantitativeCorrelationThreshold;

			FloatPair intensity=integrate(medianData.getIndices(), retentionTimes, chromatograms.get(i));

			integrationArray[i]=intensity.getOne();
			backgroundArray[i]=intensity.getTwo();
			
			// calculate trapezoidal background area
			integrationArray[i]=integrationArray[i]-backgroundArray[i];
		}

		Range range=new Range(retentionTimes[medianData.getIndices().getStart()], retentionTimes[medianData.getIndices().getStop()]);
		if (plot) {
			HashMap<String, ChartPanel> panels=new HashMap<String, ChartPanel>();
			panels.put("unnormalized", getChart(chromatograms, correlationArray, retentionTimes, range));
			panels.put("unnormalized_uncolored", getChart(chromatograms, new float[correlationArray.length], retentionTimes, range));
			panels.put("normalized", getChart(medianData.getNormalizedChromatograms(), correlationArray, retentionTimes, range));
			panels.put("median", Charter.getChart("scan", "intensity", false, toXYTrace(medianData.getMedianChromatogram(), retentionTimes, "median", null, null)));
			
			ArrayList<XYTrace> traces=getTraces(medianData.getNormalizedChromatograms(), correlationArray, retentionTimes, range);
			traces.add(0, toXYTrace(medianData.getMedianChromatogram(), retentionTimes, "median", Color.black, null, GraphType.dashedline, 6.0f));
			panels.put("traces", Charter.getChart("scan", "intensity", false, traces.toArray(new XYTrace[traces.size()])));
			
			Charter.launchCharts(peptideModSeq+" chart", panels);
		}
		
		return new TransitionRefinementData(peptideModSeq, precursorCharge, fragmentMasses, chromatograms, correlationArray, quantitativeIonsArray, integrationArray, backgroundArray, medianData.getMedianChromatogram(), range, params.getAAConstants());
	}

	public static MedianChromatogramData extractMedianChromatogramFromFixedBoundaries(float retentionTimeInSec, ArrayList<float[]> chromatograms, float[] retentionTimes, Range rtRange) {
		if (chromatograms.size()==0) return new MedianChromatogramData(new ArrayList<>(), new IntRange(0, 0), new float[0]);

		int min=retentionTimes.length-1;
		int max=0;
		for (int i = 0; i < retentionTimes.length; i++) {
			if (rtRange.contains(retentionTimes[i])) {
				min=i+1;
				break;
			}
		}
		for (int i = retentionTimes.length-1; i >=0; i--) {
			if (rtRange.contains(retentionTimes[i])) {
				max=i-1;
				break;
			}
		}
		min=Math.min(retentionTimes.length-1, min);
		max=Math.max(0, max);
		IntRange indices=new IntRange(min, max);

		// start across the entire width
		ArrayList<float[]> normalizedChromatograms=normalizeAndBackgroundSubtract(chromatograms, indices);
		
		// find the maximum point
		float[] medianChromatogram=new float[chromatograms.get(0).length];
		for (int i=0; i<medianChromatogram.length; i++) {
			TFloatArrayList list=new TFloatArrayList();
			for (float[] chromatogram : normalizedChromatograms) {
				if (chromatogram.length>i) {
					list.add(chromatogram[i]);
				} else {
					list.add(0.0f);
				}
			}
			float[] array=list.toArray();
			Arrays.sort(array);
			medianChromatogram[i]=aggregate(array);
		}
		
		MedianChromatogramData medianData=new MedianChromatogramData(normalizedChromatograms, indices, medianChromatogram);
		return medianData;
	}

	public static MedianChromatogramData extractMedianChromatogram(float retentionTimeInSec, ArrayList<float[]> chromatograms, float[] retentionTimes, Optional<float[]> maybeMedianChromatogram, boolean adjustPeakBoundaries, 
			int minNumIntegratedRTPoints, float expectedPeakWidth) {
		if (chromatograms.size()==0) return new MedianChromatogramData(new ArrayList<>(), new IntRange(0, 0), new float[0]);
		
		ArrayList<float[]> normalizedChromatograms;
		IntRange indices;
		float[] medianChromatogram;

		if (maybeMedianChromatogram.isPresent()&&maybeMedianChromatogram.get().length>0) {
			// already started with a median chromatogram
			medianChromatogram=maybeMedianChromatogram.get();
			indices=getIndexRange(retentionTimes, medianChromatogram, expectedPeakWidth, minNumIntegratedRTPoints);
			normalizedChromatograms=normalizeAndBackgroundSubtract(chromatograms, indices);
			
		} else {
			// start across the entire width
			normalizedChromatograms=normalize(chromatograms);
			
			// find the maximum point
			medianChromatogram=new float[chromatograms.get(0).length];
			for (int i=0; i<medianChromatogram.length; i++) {
				TFloatArrayList list=new TFloatArrayList();
				for (float[] chromatogram : normalizedChromatograms) {
					if (chromatogram.length>i) {
						list.add(chromatogram[i]);
					} else {
						list.add(0.0f);
					}
				}
				float[] array=list.toArray();
				Arrays.sort(array);
				medianChromatogram[i]=aggregate(array);
			}
			IntRange initialIndices=getIndexRange(retentionTimes, medianChromatogram, expectedPeakWidth, minNumIntegratedRTPoints);

			if (false) { //FIXME REMOVE TEST CODE
				String name="TEST"; //peptideModSeq+"p"+data.getPrecursorCharge();
				System.out.println("ArrayList<XYPoint> "+name+"Points=new ArrayList<>();");
				XYTrace median=new XYTrace(retentionTimes, medianChromatogram, GraphType.line, "median");
				for (XYPoint point : median.getPoints()) {
					System.out.println(name+"Points.add(new XYPoint("+point.x+", "+point.y+"));");
				}
				System.out.println("XYTrace "+name+"Trace=new XYTrace("+name+"Points, GraphType.line, \""+name+"\");");
				System.out.println("Range "+name+"Range=TransitionRefiner.getPeakRange("+name+"Trace, 25f);");
				System.out.println("XYTrace "+name+"Trimmed="+name+"Trace.trim("+name+"Range).updateColor(Color.green, 4.0f);");
				System.out.println("Charter.launchChart(\"RT\", \"Intensity\", true, "+name+"Trace, "+name+"Trimmed);");
			}
			
			if (!adjustPeakBoundaries) {
				
				Range testRange=new Range(retentionTimes[initialIndices.getStart()], retentionTimes[initialIndices.getStop()]);
				if (!testRange.contains(retentionTimeInSec)) {
						// if it wasn't inferred and our boundaries were outside the range, we need to reset the range 
					int index=Arrays.binarySearch(retentionTimes, retentionTimeInSec);

					if (index < 0) {
						index = -(index + 1);
					}

					// Issue #79. Guard against out-of-bounds. This is NOT an `else`
					// clause of the `if` above, because the insertion point returned
					// by binary search can be the end of the array.
					// Note that the resulting index could still be out of bounds, but
					// only if the array is length zero. However, it's probably safe to
					// assume that the arrays have matching length, and the call to
					// `getIndexRange()` above would have failed if the length was zero.
					if (index >= retentionTimes.length) {
						index = retentionTimes.length - 1;
					}

					initialIndices=getIndexRange(retentionTimes, medianChromatogram, index, expectedPeakWidth, minNumIntegratedRTPoints);
				}
				testRange=new Range(retentionTimes[initialIndices.getStart()], retentionTimes[initialIndices.getStop()]);
			}

			// then refine on the local area
			normalizedChromatograms=normalizeAndBackgroundSubtract(chromatograms, initialIndices);
			
			// find the maximum point
			medianChromatogram=new float[chromatograms.get(0).length];
			for (int i=0; i<medianChromatogram.length; i++) {
				TFloatArrayList list=new TFloatArrayList();
				for (float[] chromatogram : normalizedChromatograms) {
					if (chromatogram.length>i) {
						list.add(chromatogram[i]);
					} else {
						list.add(0.0f);
					}
				}
				float[] array=list.toArray();
				Arrays.sort(array);
				medianChromatogram[i]=aggregate(array);
			}
			indices=getIndexRange(retentionTimes, medianChromatogram, expectedPeakWidth, minNumIntegratedRTPoints);
		}

		MedianChromatogramData medianData=new MedianChromatogramData(normalizedChromatograms, indices, medianChromatogram);
		return medianData;
	}
	
	private static float aggregate(float[] values) {
		if (false) {
			return General.mean(values);
		} else {
			return QuickMedian.median(values);
		}
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
			for (int ionIndex=0; ionIndex<allChromatograms.length; ionIndex++) {
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

	public static float calculateCorrelation(IntRange indices, float[] medianChromatogram, float[] normalizedChromatogram) {
		float medianMean=General.mean(medianChromatogram, indices.getStart(), indices.getStop());
		return calculateCorrelation(medianMean, indices, medianChromatogram, normalizedChromatogram);
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
			if (denominator==0.0f) return Float.MIN_VALUE;
			
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

	public static Range getPeakRange(XYTrace trace, int minNumIntegratedRTPoints, float expectedPeakWidth) {
		Pair<double[], double[]> pair=trace.toArrays();
		float[] x=General.toFloatArray(pair.x);
		float[] y=General.toFloatArray(pair.y);
		IntRange indicies= getIndexRange(x, y, expectedPeakWidth, minNumIntegratedRTPoints);
		return new Range(trace.getPoints().get(indicies.getStart()).x, trace.getPoints().get(indicies.getStop()).x);
	}
	
	public static IntRange getIndexRange(float[] x, float[] y, float expectedPeakWidth, int minNumIntegratedRTPoints) {
		// default to be the middle of the window
		int maxIndex=(int)Math.floor((y.length-1)/2.0f);
		if (maxIndex<0) maxIndex=0;
		float maxY=y[maxIndex];
		
		for (int i = 0; i < y.length; i++) {
			if (y[i]>maxY) {
				maxIndex=i;
				maxY=y[i];
			}
		}
		return getIndexRange(x, y, maxIndex, expectedPeakWidth, minNumIntegratedRTPoints);
	}
	
	public static IntRange getIndexRange(float[] x, float[] y, int maxIndex, float expectedPeakWidth, int minNumIntegratedRTPoints) {
		float fiftyPercentPoint=y[maxIndex]/2.0f;
		float threshold=y[maxIndex]*MINIMUM_THRESHOLD_PERCENTAGE; // 1% of max
		float rtThreshold=expectedPeakWidth*PERCENT_RT_UPHILL_PEAKWIDTH; // if we start moving by more than 20% of the expected peakwidth
		float altRTThreshold=2.1f*QuickMedian.median(General.firstDerivative(x));
		if (altRTThreshold>rtThreshold) {
			// allow at least 2 points of deviation!
			rtThreshold=altRTThreshold;
		}
		
		// left of center (decreasing index)
		int firstData=maxIndex;
		int lowestIndex=maxIndex;
		for (int i=maxIndex-1; i>=0; i--) {
			// navigate down the slope of the peak:
			if (y[i]<y[lowestIndex]) {
				lowestIndex=i;
			}
			
			// if we're more than 10% from the lowest point and we're less than 50% of the max
			//System.out.println("TEST LEFT 50% "+i+"("+x[i]+") --> "+y[i]+" ("+Math.abs(x[i]-x[lowestIndex])+" > "+rtThreshold+") "+(Math.abs(x[i]-x[lowestIndex])>rtThreshold&&fiftyPercentPoint>y[lowestIndex]));
			if (Math.abs(x[i]-x[lowestIndex])>rtThreshold&&fiftyPercentPoint>y[lowestIndex]) {
				firstData=lowestIndex; // switch back to actual lowest point
				//System.out.println("HIT LEFT 50% "+i+"("+x[i]+") --> "+y[i]+" ("+Math.abs(x[i]-x[lowestIndex])+" > "+rtThreshold+") "+(Math.abs(x[i]-x[lowestIndex])>rtThreshold&&fiftyPercentPoint>y[lowestIndex]));
				break;
			}
			
			// if we're lower than the previous local minimum, set the new minimum
			if (y[i]<y[firstData]) {
				firstData=i;
			}
			
			// create peak boundary if the local minimum is less than 1%
			if (y[firstData]<threshold) {
				firstData=i-1; // go out one more spot
				//System.out.println("HIT LEFT THRESHOLD "+firstData+" --> "+y[firstData]);
				break;
			}
		}
		
		// right of center (increasing index)
		lowestIndex=maxIndex;
		int lastData=maxIndex;
		for (int i=maxIndex+1; i<y.length; i++) {
			// navigate down the slope of the peak:
			// count the number of consecutive uphill points (count the aggregate, so +1 for increasing, -1 for decreasing)
			if (y[i]<y[lowestIndex]) {
				lowestIndex=i;
			}
			
			// if we're more than 10% from the lowest point and we're less than 50% of the max
			//System.out.println("TEST RIGHT 50% "+i+"("+x[i]+") --> "+y[i]+" ("+Math.abs(x[i]-x[lowestIndex])+" > "+rtThreshold+") "+(Math.abs(x[i]-x[lowestIndex])>rtThreshold&&fiftyPercentPoint>y[lowestIndex]));
			
			if (Math.abs(x[i]-x[lowestIndex])>rtThreshold&&fiftyPercentPoint>y[lowestIndex]) {
				lastData=lowestIndex; // switch back to actual lowest point
				//System.out.println("HIT RIGHT 50% "+i+"("+x[i]+") --> "+y[i]+" ("+Math.abs(x[i]-x[lowestIndex])+" > "+rtThreshold+") "+(Math.abs(x[i]-x[lowestIndex])>rtThreshold&&fiftyPercentPoint>y[lowestIndex]));
				break;
			}
			
			// if we're lower than the previous local minimum, set the new minimum
			if (y[i]<y[lastData]) {
				lastData=i;
			}
			
			// create peak boundary if the local minimum is less than 1%
			if (y[lastData]<threshold) {
				lastData=i+1; // go out one more spot
				//System.out.println("HIT RIGHT THRESHOLD "+lastData+" --> "+y[lastData]);
				break;
			}
		}

		
		int startIndex=(firstData<=0)?(0):(firstData);
		int stopIndex=(lastData>=y.length-1)?(y.length-1):(lastData);
		
		boolean lastUpdatedStop=true;
		// increment by checking for the best value on either side
		while (stopIndex-startIndex<(minNumIntegratedRTPoints-1)) {
			float startSide=NO_AVAILABLE_INTENSITY_MARKER;
			float stopSide=NO_AVAILABLE_INTENSITY_MARKER;
			
			if (startIndex>0) {
				startSide=y[startIndex-1];
			}
			if (stopIndex<y.length-1) {
				stopSide=y[stopIndex+1];
			}

			if (startSide>stopSide) {
				startIndex=startIndex-1;
				lastUpdatedStop=false;
			} else if (startSide<stopSide) {
				stopIndex=stopIndex+1;
			} else if (startSide==NO_AVAILABLE_INTENSITY_MARKER) {
				// both startSide and stopSide==NO_AVAILABLE_INTENSITY_MARKER
				break;
			} else {
				// if there's no difference, then choose whatever was not last chosen
				if (lastUpdatedStop) {
					startIndex=startIndex-1;
					lastUpdatedStop=false;
				} else {
					stopIndex=stopIndex+1;
					lastUpdatedStop=true;
				}
			}
		}
		IntRange indices=new IntRange(startIndex, stopIndex);
		return indices;
	}
	
	public static HashMap<String, ChartPanel> getChartPanels(TransitionRefinementData data) {
		HashMap<String, ChartPanel> panels=new HashMap<String, ChartPanel>();
		float[] rts=null;
		if (data!=null&&data.getRtArray()!=null&&data.getRtArray().isPresent()) {
			rts=data.getRtArray().get();
		}
		XYTrace median = toXYTrace(data.getMedianChromatogram(), rts, "median", null, data.getRange());
		
		float max=General.max(data.getMedianChromatogram());
		XYTrace start=new XYTrace(new float[] {data.getRange().getStart()/60f, data.getRange().getStart()/60f}, new float[] {0, max},GraphType.line, "leftBoundary", Color.black, null);
		XYTrace stop=new XYTrace(new float[] {data.getRange().getStop()/60f, data.getRange().getStop()/60f}, new float[] {0, max},GraphType.line, "rightBoundary", Color.black, null);

		panels.put("unnormalized", getChart(data.getChromatograms(), data.getCorrelationArray(), rts, data.getRange()));
		panels.put("median", Charter.getChart("Retention Time (min)", "intensity", false, median, start, stop));
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
