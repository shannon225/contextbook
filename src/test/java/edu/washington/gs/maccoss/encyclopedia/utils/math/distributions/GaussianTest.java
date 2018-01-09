package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;
import junit.framework.TestCase;

public class GaussianTest extends TestCase {
	// calculates error in trapezoidal areas given N number of points across the peak 
	public static void main(String[] args) {
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
				float x=r.getStart()+offset*increment+interval;
				
				float totalArea=0.0f;
				for (int i=1; i<numPointsAcrossPeak; i++) {
					if (!r.contains(x)) continue;
					float y=(float)dist.getProbability(x);
					
					float trapezoidalArea=(x-lastX)*((lastY+y)/2.0f);
					totalArea+=trapezoidalArea;
					lastX=x;
					lastY=y;
					x=x+interval;
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
