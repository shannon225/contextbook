package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Collections;

import gnu.trove.map.hash.TDoubleObjectHashMap;
import junit.framework.TestCase;

public class RangeTest extends TestCase {
	public void testInterpolation() {
		Range r=new Range(30f, 40f);
		// middle
		assertEquals(500f, r.linearInterp(35f, 0f, 1000f));
		assertEquals(45f, r.linearInterp(35f, 40f, 50f));
		
		// boundaries
		assertEquals(40f, r.linearInterp(30f, 40f, 50f));
		assertEquals(50f, r.linearInterp(40f, 40f, 50f));
		

		// middle
		assertEquals(35f, r.mapBackToRange(500f, 0f, 1000f));
		assertEquals(35f, r.mapBackToRange(45f, 40f, 50f));
		
		// boundaries
		assertEquals(30f, r.mapBackToRange(40f, 40f, 50f));
		assertEquals(40f, r.mapBackToRange(50f, 40f, 50f));
	}

	public void testWidestRange() {
		ArrayList<Range> ranges=new ArrayList<Range>();
		ranges.add(new Range(10f, 20f));
		ranges.add(new Range(30f, 40f));
		ranges.add(new Range(-10f, 20f));
		ranges.add(new Range(10f, 20f));
		
		Range widest=Range.getWidestRange(ranges);
		assertEquals(-10f, widest.getStart());
		assertEquals(40f, widest.getStop());
	}
	
	public void testContains() {
		Range r1=new Range(100, 200);
		Range r2=new Range(199, 199);
		Range r3=new Range(199, 201);
		Range r4=new Range(200, 300);
		Range r5=new Range(300, 400);

		Range r6=new Range(101, 101);
		Range r7=new Range(99, 101);
		Range r8=new Range(98, 99);
		
		Range r9=new Range(0, 300);
		
		assertTrue(r1.contains(r2));
		assertTrue(r1.contains(r3));
		assertTrue(r1.contains(r4));
		assertFalse(r1.contains(r5));
		assertTrue(r1.contains(r6));
		assertTrue(r1.contains(r7));
		assertFalse(r1.contains(r8));
		assertTrue(r1.contains(r9));
	}
	
	public void testComparator() {
		TDoubleObjectHashMap<Range> precursorIsolationRanges=new TDoubleObjectHashMap<>();
	    precursorIsolationRanges.put(415, new Range(400,430));
	    precursorIsolationRanges.put(439, new Range(429,449));
	    precursorIsolationRanges.put(458, new Range(448,468));
	    precursorIsolationRanges.put(477, new Range(467,487));
	    precursorIsolationRanges.put(491, new Range(486,496));
	    precursorIsolationRanges.put(510, new Range(505,515));
	    precursorIsolationRanges.put(519, new Range(514,524));
	    precursorIsolationRanges.put(528, new Range(523,533));
	    precursorIsolationRanges.put(537, new Range(532,542));
	    precursorIsolationRanges.put(546, new Range(541,551));
	    precursorIsolationRanges.put(555, new Range(550,560));
	    precursorIsolationRanges.put(564, new Range(559,569));
	    precursorIsolationRanges.put(573, new Range(568,578));
	    precursorIsolationRanges.put(582, new Range(577,587));
	    precursorIsolationRanges.put(591, new Range(586,596));
	    precursorIsolationRanges.put(600, new Range(595,605));
	    precursorIsolationRanges.put(609, new Range(604,614));
	    precursorIsolationRanges.put(618, new Range(613,623));
	    precursorIsolationRanges.put(627, new Range(622,632));
	    precursorIsolationRanges.put(636, new Range(631,641));
	    precursorIsolationRanges.put(645, new Range(640,650));
	    precursorIsolationRanges.put(654, new Range(649,659));
	    precursorIsolationRanges.put(663, new Range(658,668));
	    precursorIsolationRanges.put(672, new Range(667,677));
	    precursorIsolationRanges.put(681, new Range(676,686));
	    precursorIsolationRanges.put(690, new Range(685,695));
	    precursorIsolationRanges.put(699, new Range(694,704));
	    precursorIsolationRanges.put(708, new Range(703,713));
	    precursorIsolationRanges.put(717, new Range(712,722));
	    precursorIsolationRanges.put(726, new Range(721,731));
	    precursorIsolationRanges.put(735, new Range(730,740));
	    precursorIsolationRanges.put(744, new Range(739,749));
	    precursorIsolationRanges.put(755.5, new Range(748,763));
	    precursorIsolationRanges.put(764.5, new Range(757,772));
	    precursorIsolationRanges.put(781, new Range(771,791));
	    precursorIsolationRanges.put(795, new Range(785,805));
	    precursorIsolationRanges.put(814, new Range(804,824));
	    precursorIsolationRanges.put(833, new Range(823,843));
	    precursorIsolationRanges.put(852, new Range(842,862));
	    precursorIsolationRanges.put(871, new Range(861,881));
	    precursorIsolationRanges.put(890, new Range(880,900));
	    precursorIsolationRanges.put(911.5, new Range(899,924));
	    precursorIsolationRanges.put(930.5, new Range(918,943));
	    precursorIsolationRanges.put(957, new Range(942,972));
	    precursorIsolationRanges.put(981, new Range(966,996));
	    precursorIsolationRanges.put(1015, new Range(995,1035));
	    precursorIsolationRanges.put(1044, new Range(1024,1064));
	    precursorIsolationRanges.put(1083, new Range(1063,1103));
	    precursorIsolationRanges.put(1124.5, new Range(1102,1147));
	    precursorIsolationRanges.put(1163.5, new Range(1141,1186));
	    precursorIsolationRanges.put(1207.5, new Range(1185,1230));
	    precursorIsolationRanges.put(1254, new Range(1229,1279));
	    precursorIsolationRanges.put(1298, new Range(1273,1323));
	    precursorIsolationRanges.put(1347, new Range(1322,1372));
	    precursorIsolationRanges.put(1396, new Range(1371,1421));
	    precursorIsolationRanges.put(1445, new Range(1420,1470));
	    precursorIsolationRanges.put(1494, new Range(1469,1519));
	    precursorIsolationRanges.put(1543, new Range(1518,1568));
	    precursorIsolationRanges.put(1594.5, new Range(1567,1622));
	    precursorIsolationRanges.put(1645.5, new Range(1616,1675));
	    
	    ArrayList<Range> ranges=new ArrayList<>(precursorIsolationRanges.valueCollection());
	    Collections.sort(ranges);
	    
	    for (double center : precursorIsolationRanges.keys()) {
	    		double target=center+Math.random()/10.0; // shift by <0.1 to make targets different than centers without extending beyond expected boundaries
	    		Range key=new Range(target, target);
	    		int nearestIndex=Collections.binarySearch(ranges, key, Range.RANGE_CONTAINS_COMPARATOR);
	    		assert(nearestIndex>=0);
	    		
			assertEquals(precursorIsolationRanges.get(center), ranges.get(nearestIndex));
		}
	}
}
