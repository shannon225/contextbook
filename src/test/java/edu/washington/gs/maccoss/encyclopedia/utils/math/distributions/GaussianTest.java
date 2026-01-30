package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearInterpolatedFunction;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LogQuadraticInterpolatedFunction;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import junit.framework.TestCase;

public class GaussianTest extends TestCase {
	// calculates error in trapezoidal areas given N number of points across the peak 
	public static void main2(String[] args) {
		Gaussian dist=new Gaussian(0, 1, 1);
		
		Range r=new Range(-3, 3);
		float baseline=(float)(dist.getCDF(r.getStop())-dist.getCDF(r.getStart()));
		for (int numPointsAcrossPeak=35; numPointsAcrossPeak>0; numPointsAcrossPeak--) {
			float interval=r.getRange()/(float)numPointsAcrossPeak;
			float increment=interval/100f;

			TFloatArrayList totalAreas=new TFloatArrayList();
			for (int offset=0; offset<100; offset++) {
				float lastX=r.getStart()+offset*increment;
				float lastY=0.0f;
				float xi=r.getStart()+offset*increment+interval;
				
				TFloatArrayList xs=new TFloatArrayList();
				TFloatArrayList ys=new TFloatArrayList();
				for (int i=0; i<numPointsAcrossPeak; i++) {
					float y=(float)dist.getProbability(xi);
					//if (Math.random()<0.1) y=0.0f;
					xs.add(xi);
					ys.add(y);
					xi=xi+interval;
				}
				//ys=new TFloatArrayList(SkylineSGFilter.paddedSavitzkyGolaySmooth(ys.toArray()));
				
				float x=0;
				float totalArea=0.0f;
				for (int i=0; i<xs.size(); i++) {
					x=xs.get(i);
					float y=ys.get(i);
					if (!r.contains(x)) continue;
					
					float trapezoidalArea=(x-lastX)*((lastY+y)/2.0f);
					totalArea+=trapezoidalArea;
					lastX=x;
					lastY=y;
				}
				float trapizoidalArea=(x-lastX)*((lastY)/2.0f);
				totalArea+=trapizoidalArea;
				totalAreas.add(totalArea-baseline);
			}
			
			final float[] areas=totalAreas.toArray();
			Arrays.sort(areas);
			final float mean=General.mean(areas);
			System.out.println(numPointsAcrossPeak+"\t"+100*mean+"\t"+100*areas[Math.round(areas.length*0.05f)]+"\t"+100*areas[Math.round(areas.length*0.95f)]);
		}
	}
	
	public static void main3(String[] args) {
		Gaussian dist=new Gaussian(0, 1, 1);
		Range r=new Range(-3, 3);

		double noiseLevel=0.00;
		double noisePercentage=0.0;
		int seed=RandomGenerator.randomInt(95);
		
		int nPointsAcrossPeak=2; // 5 with zeros
		
		double peakSpacing=r.getRange()/nPointsAcrossPeak;

		seed=RandomGenerator.randomInt(seed);
		double offset=RandomGenerator.random(seed)*peakSpacing;

		TDoubleArrayList xList=new TDoubleArrayList(nPointsAcrossPeak+2);
		TDoubleArrayList yList=new TDoubleArrayList(nPointsAcrossPeak+2);

		
		// add the peaks in between
		for (int i=0; i<nPointsAcrossPeak; i++) {
			double x=r.getStart()+offset+i*peakSpacing;
			double y=dist.getProbability(x);
			
			double noisePercent=noiseLevel*RandomGenerator.random(seed+1)-noiseLevel/2.0f;
			if (RandomGenerator.random(seed+2)<noisePercentage) {
				y=Math.max(0.0f, y+noisePercent*y);
			}
			
			xList.add(x);
			yList.add(y);
		}
		// add a peak before that's 0
		xList.insert(0, r.getStart()+offset-peakSpacing);
		yList.insert(0, 0.000);
		// add a peak after that's 0
		xList.add(r.getStart()+offset+(nPointsAcrossPeak)*peakSpacing);
		yList.add(0.000);

		double[] xs=xList.toArray();
		double[] ys=yList.toArray();
		
		LinearInterpolatedFunction linear=new LinearInterpolatedFunction(xs, ys);
		LogQuadraticInterpolatedFunction quad=new LogQuadraticInterpolatedFunction(xs, ys);
		
		System.out.println(General.min(xs)+", "+General.max(xs));
		for (float i=-4.0f; i<+4.0f; i=i+0.1f) {
			System.out.println(i+"\t"+dist.getProbability(i)+"\t"+linear.getYValue(i)+"\t"+quad.getYValue(i));
		}
		
		System.out.println("\nIntegrations lin:"+linear.integrate(r.getStart(), r.getStop())+", quad:"+quad.integrate(r.getStart(), r.getStop()));
		for (int i=0; i<xs.length; i++) {
			System.out.println(xs[i]+"\t"+ys[i]);
		}
	}
	
	// introducing error data
	public static void mainNoise(String[] args) {
		boolean isLinear=false;
		boolean isApex=false;
		boolean isRT=false;
		
		double mean=0.0;
		double stdev=1.0;
		Gaussian dist=new Gaussian(mean, stdev, 1.0);
		Range r=new Range(-3, 3);

		int nTrials=100000;
		int seed=RandomGenerator.randomInt(2);
		
		float baseline=(float)(dist.getCDF(r.getStop())-dist.getCDF(r.getStart()));
		double apex=dist.getPDF(dist.getMean());
		
		for (float noiseLevel=0; noiseLevel<=0.5; noiseLevel+=0.01f) {
		
			int minPeaks=4; // 6 with 0s
			int maxPeaks=8; // 12 with 0s
			
			TFloatArrayList[] allErrors=new TFloatArrayList[maxPeaks-minPeaks+1];
			for (int i=0; i<allErrors.length; i++) {
				allErrors[i]=new TFloatArrayList();
			}
			for (int nPointsAcrossPeak=minPeaks; nPointsAcrossPeak<=maxPeaks; nPointsAcrossPeak++) {
				for (int trials=0; trials<nTrials; trials++) {
					float peakSpacing=r.getRange()/nPointsAcrossPeak;
	
					seed=RandomGenerator.randomInt(seed);
					float offset=RandomGenerator.random(seed)*peakSpacing;
	
					TFloatArrayList xList=new TFloatArrayList(nPointsAcrossPeak+2);
					TFloatArrayList yList=new TFloatArrayList(nPointsAcrossPeak+2);
					
					// add the peaks in between
					int actualPointsAbove1p=0;
					for (int i=0; i<nPointsAcrossPeak; i++) {
						seed=RandomGenerator.randomInt(seed);
						
						float x=r.getStart()+offset+i*peakSpacing;
						float y=(float)dist.getProbability(x);
						
						// +/- up to 100% error
						float noisePercent=2*noiseLevel*RandomGenerator.random(seed+nPointsAcrossPeak)-noiseLevel;
						y=Math.max(0.0f, y+noisePercent*(float)y);
						
						if (y>0.01f*apex) {
							actualPointsAbove1p++;
						}
						
						xList.add(x);
						yList.add(y);
					}
					// add a peak before that's 0
					xList.insert(0, r.getStart()+offset-peakSpacing);
					yList.insert(0, Math.min(0.00f, yList.get(0)));
					// add a peak after that's 0
					xList.add(r.getStart()+offset+(nPointsAcrossPeak)*peakSpacing);
					yList.add(Math.min(0.00f, yList.get(yList.size()-1)));
	
					float[] xs=xList.toArray();
					float[] ys=yList.toArray();
					
					float percentError;
					if (isApex) {
						if (isRT) {
							float observedRT;
							if (isLinear) {
								LinearInterpolatedFunction linear=new LinearInterpolatedFunction(General.toDoubleArray(xs), General.toDoubleArray(ys));
								observedRT=(float)linear.getApex(r.getStart(), r.getStop()).x;
							} else {
								LogQuadraticInterpolatedFunction quad=new LogQuadraticInterpolatedFunction(General.toDoubleArray(xs), General.toDoubleArray(ys));
								observedRT=(float)quad.getApex(r.getStart(), r.getStop()).x;
							}
							
							float deltaRT=(float)Math.abs(observedRT-mean);
							
							percentError=deltaRT;
						} else {
							float totalArea;
							if (isLinear) {
								LinearInterpolatedFunction linear=new LinearInterpolatedFunction(General.toDoubleArray(xs), General.toDoubleArray(ys));
								totalArea=(float)linear.getApex(r.getStart(), r.getStop()).y;
							} else {
								LogQuadraticInterpolatedFunction quad=new LogQuadraticInterpolatedFunction(General.toDoubleArray(xs), General.toDoubleArray(ys));
								totalArea=(float)quad.getApex(r.getStart(), r.getStop()).y;
							}
							
							percentError=100.0f*totalArea/(float)apex;
						}
					} else {
						float totalArea;
						if (isLinear) {
							LinearInterpolatedFunction linear=new LinearInterpolatedFunction(General.toDoubleArray(xs), General.toDoubleArray(ys));
							totalArea=(float)linear.integrate(r.getStart(), r.getStop());
						} else {
							LogQuadraticInterpolatedFunction quad=new LogQuadraticInterpolatedFunction(General.toDoubleArray(xs), General.toDoubleArray(ys));
							totalArea=(float)quad.integrate(r.getStart(), r.getStop());
						}
						
						percentError=100.0f*totalArea/baseline;
					}
					
					if (actualPointsAbove1p>=minPeaks&&actualPointsAbove1p<=maxPeaks) {
						allErrors[actualPointsAbove1p-minPeaks].add(Float.valueOf(percentError));
					}
				}
			}
			
			Float[][] allErrorArray=new Float[allErrors.length][];

			System.out.print(noiseLevel);
			for (int i=0; i<allErrorArray.length; i++) {
				float median=QuickMedian.median(allErrors[i].toArray());
				System.out.print("\t");
				System.out.print(median);
			}
			System.out.println();
		}
	}
	

	// boxplot data
	public static void main(String[] args) {
		boolean isLinear=true;
		boolean isApex=true;
		boolean isRT=true;
		
		double mean=0.0;
		double stdev=1.0;
		Gaussian dist=new Gaussian(mean, stdev, 1.0);
		Range r=new Range(-3, 3);

		float noiseLevel=0.0f;
		float noisePercentage=0f;
		int nTrials=100000;
		int seed=RandomGenerator.randomInt(2);
		
		float baseline=(float)(dist.getCDF(r.getStop())-dist.getCDF(r.getStart()));
		double apex=dist.getPDF(dist.getMean());
		
		int minPeaks=4; // 6 with 0s
		int maxPeaks=10; // 12 with 0s
		
		ArrayList[] allErrors=new ArrayList[maxPeaks-minPeaks+1];
		for (int i=0; i<allErrors.length; i++) {
			allErrors[i]=new ArrayList<Float>();
		}
		for (int nPointsAcrossPeak=minPeaks; nPointsAcrossPeak<=maxPeaks; nPointsAcrossPeak++) {
			for (int trials=0; trials<nTrials; trials++) {
				float peakSpacing=r.getRange()/nPointsAcrossPeak;

				seed=RandomGenerator.randomInt(seed);
				float offset=RandomGenerator.random(seed)*peakSpacing;

				TFloatArrayList xList=new TFloatArrayList(nPointsAcrossPeak+2);
				TFloatArrayList yList=new TFloatArrayList(nPointsAcrossPeak+2);
				
				// add the peaks in between
				int actualPointsAbove1p=0;
				for (int i=0; i<nPointsAcrossPeak; i++) {
					float x=r.getStart()+offset+i*peakSpacing;
					float y=(float)dist.getProbability(x);
					
					float noisePercent=noiseLevel*RandomGenerator.random(seed+1)-noiseLevel/2.0f;
					if (RandomGenerator.random(seed+2)<noisePercentage) {
						y=Math.max(0.0f, y+noisePercent*y);
					}
					if (y>0.01f*0.39894228f) {
						actualPointsAbove1p++;
					}
					
					xList.add(x);
					yList.add(y);
				}
				// add a peak before that's 0
				xList.insert(0, r.getStart()+offset-peakSpacing);
				yList.insert(0, Math.min(0.00f, yList.get(0)));
				// add a peak after that's 0
				xList.add(r.getStart()+offset+(nPointsAcrossPeak)*peakSpacing);
				yList.add(Math.min(0.00f, yList.get(yList.size()-1)));

				float[] xs=xList.toArray();
				float[] ys=yList.toArray();
				
				float percentError;
				if (isApex) {
					if (isRT) {
						float observedRT;
						if (isLinear) {
							LinearInterpolatedFunction linear=new LinearInterpolatedFunction(General.toDoubleArray(xs), General.toDoubleArray(ys));
							observedRT=(float)linear.getApex(r.getStart(), r.getStop()).x;
						} else {
							LogQuadraticInterpolatedFunction quad=new LogQuadraticInterpolatedFunction(General.toDoubleArray(xs), General.toDoubleArray(ys));
							observedRT=(float)quad.getApex(r.getStart(), r.getStop()).x;
						}
						
						float deltaRT=(float)Math.abs(observedRT-mean);
						
						percentError=deltaRT;
					} else {
						float totalArea;
						if (isLinear) {
							LinearInterpolatedFunction linear=new LinearInterpolatedFunction(General.toDoubleArray(xs), General.toDoubleArray(ys));
							totalArea=(float)linear.getApex(r.getStart(), r.getStop()).y;
						} else {
							LogQuadraticInterpolatedFunction quad=new LogQuadraticInterpolatedFunction(General.toDoubleArray(xs), General.toDoubleArray(ys));
							totalArea=(float)quad.getApex(r.getStart(), r.getStop()).y;
						}
						
						percentError=100.0f*totalArea/(float)apex;
					}
				} else {
					float totalArea;
					if (isLinear) {
						LinearInterpolatedFunction linear=new LinearInterpolatedFunction(General.toDoubleArray(xs), General.toDoubleArray(ys));
						totalArea=(float)linear.integrate(r.getStart(), r.getStop());
					} else {
						LogQuadraticInterpolatedFunction quad=new LogQuadraticInterpolatedFunction(General.toDoubleArray(xs), General.toDoubleArray(ys));
						totalArea=(float)quad.integrate(r.getStart(), r.getStop());
					}
					
					percentError=100.0f*totalArea/baseline;
				}
				
				if (actualPointsAbove1p>=minPeaks&&actualPointsAbove1p<=maxPeaks) {
					allErrors[actualPointsAbove1p-minPeaks].add(Float.valueOf(percentError));
				}
			}
		}
		
		Float[][] allErrorArray=new Float[allErrors.length][];
		for (int i=0; i<allErrorArray.length; i++) {
			allErrorArray[i]=(Float[])allErrors[i].toArray(new Float[0]);
		}
		
		allErrorArray=General.transposeMatrix(allErrorArray);

		for (int nPointsAcrossPeak=minPeaks; nPointsAcrossPeak<=maxPeaks; nPointsAcrossPeak++) {
			if (nPointsAcrossPeak>minPeaks) {
				System.out.print("\t");
			}
			System.out.print(nPointsAcrossPeak);
		}
		System.out.println();
		for (int i=0; i<allErrorArray.length; i++) {
			for (int j=0; j<allErrorArray[i].length; j++) {
				if (j>0) {
					System.out.print("\t");
				}
				if (allErrorArray[i][j]!=null) {
					System.out.print(allErrorArray[i][j]);
				}
			}
			System.out.println();
		}
	}
	
	public void testGaussian() {
		float prior=7f;
		Distribution g=new Gaussian(0, 1, prior);
		float[] xs=new float[] {-3, -2, -1, 0, 1, 2, 3};
		float[] cdfs=new float[] {0.0013499672813147567f, 0.022750062887256395f, 0.15865526383236372f, 0.500000000f, 0.8413447361676363f, 0.9772499371127437f, 0.9986500327186852f};
		float[] pdfs=new float[] {0.0044318484119380075f, 0.05399096651318806f, 0.24197072451914337f, 0.3989422804014327f, 0.24197072451914337f, 0.05399096651318806f, 0.0044318484119380075f};
		for (int i=0; i<xs.length; i++) {
			assertEquals(cdfs[i], g.getCDF(xs[i]), 0.00001f);
			assertEquals(pdfs[i], g.getPDF(xs[i]), 0.00001f);
			assertEquals(pdfs[i]*prior, g.getProbability(xs[i]), 0.00001f);
		}
	}

}
