package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import edu.washington.gs.maccoss.encyclopedia.utils.SparseIndexMap;
import junit.framework.TestCase;

public class SparseXCorrSpectrumTest extends TestCase {
	public void testDotProduct() {
		SparseIndexMap m1=new SparseIndexMap();
		m1.adjustOrPutValue(100, 5, 5);
		m1.adjustOrPutValue(200, 4, 4);
		m1.adjustOrPutValue(1000, 5, 5); // stacked
		m1.adjustOrPutValue(1000, 6, 6);
		m1.adjustOrPutValue(10000, 8, 8);
		m1.adjustOrPutValue(99999900, 9, 9);
		m1.adjustOrPutValue(99999901, 7, 7);
		SparseXCorrSpectrum s1=new SparseXCorrSpectrum(m1, 2000.0, 0.01f, 100000000);

		SparseIndexMap m2=new SparseIndexMap();
		m2.adjustOrPutValue(110, 5, 5);
		m2.adjustOrPutValue(240, 4, 4);
		m2.adjustOrPutValue(1500, 5, 5);
		m2.adjustOrPutValue(1000, 6, 6);
		m2.adjustOrPutValue(10000, 8, 8);
		m2.adjustOrPutValue(99999900, 9, 9);
		m2.adjustOrPutValue(99499901, 7, 7);
		SparseXCorrSpectrum s2=new SparseXCorrSpectrum(m2, 2000.0, 0.01f, 100000000);

		assertEquals(11*6+8*8+9*9, s1.dotProduct(s2), 0.001f);
	}
	
	public void testToArray() {
		SparseIndexMap m1=new SparseIndexMap();
		m1.adjustOrPutValue(100, 5, 5);
		m1.adjustOrPutValue(200, 4, 4);
		m1.adjustOrPutValue(1000, 5, 5); // stacked
		m1.adjustOrPutValue(1000, 6, 6);
		m1.adjustOrPutValue(10000, 8, 8);
		m1.adjustOrPutValue(9999900, 9, 9);
		m1.adjustOrPutValue(9999901, 7, 7);
		SparseXCorrSpectrum s1=new SparseXCorrSpectrum(m1, 2000.0, 0.01f, 10000000);
		
		float[] array=s1.toArray();
		
		assertEquals(10000000, array.length);
		for (int i=0; i<array.length; i++) {
			if (m1.contains(i)) {
				assertEquals(m1.getPeak(i).intensity, array[i], 0.001f);
			} else {
				assertEquals(0.0f, array[i], 0.001f);
			}
		}
	}
}
