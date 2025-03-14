package edu.washington.gs.maccoss.encyclopedia.utils.graphing;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import junit.framework.TestCase;

public class XYTraceTest extends TestCase {
	public void testRound() {
		double[] x = new double[] { 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.49, 2.7, 2.9 };
		double[] y = new double[] { 11, 12, 13, 14, 15, 6, 17, 18, 119, 110, 111, 112, 113, 114, 115 };

		XYTrace trace = new XYTrace(x, y, GraphType.area, "name");
		XYTrace round = XYTrace.round(trace, 1);

		Pair<double[], double[]> xy=round.toArrays();
		float[] expectedx = new float[] { 0.0f, 1.0f, 2.0f, 3.0f };
		float[] expectedy = new float[] { 11.5f, 13.0f, 97.16666666666667f, 114.5f };

		assertEquals(expectedx.length, xy.x.length);
		assertEquals(expectedy.length, xy.y.length);
		for (int i = 0; i < expectedy.length; i++) {
			assertEquals(expectedx[i], xy.x[i], 0.0001);
			assertEquals(expectedy[i], xy.y[i], 0.0001);
		}
	}
	
	public void testIntegrate() {
		double[] x = new double[] { 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.49, 2.7, 2.9 };
		double[] y = new double[] { 11, 12, 13, 14, 15, 6, 17, 18, 119, 110, 111, 112, 113, 114, 115 };

		XYTrace trace = new XYTrace(x, y, GraphType.area, "name");
		

		assertEquals(101.485, trace.integrate(new Range(2, 2.9)), 0.0001);
		assertEquals(101.485, trace.integrate(new Range(2, 100000)), 0.0001);
		assertEquals(114.5*0.2, trace.integrate(new Range(2.7, 2.9)), 0.0001);
		//assertEquals(11.5*.2, trace.integrate(new Range(0.0, 0.3)), 0.0001);
		//assertEquals(11.5*.2, trace.integrate(new Range(0.1, 0.3)), 0.0001);
		//assertEquals(0.0, trace.integrate(new Range(12.0, 15.0)));
		//assertEquals(0.0, trace.integrate(new Range(-10, 0.1)));
	}
	public void testCorrelation() {

		double[] ax = new double[] { 1, 3, 5, 7, 9 };
		double[] ay = new double[] { 0, 5, 10, 5, 0 };

		XYTrace traceA = new XYTrace(ax, ay, GraphType.area, "name");


		double[] bx = new double[] { 2, 4, 6, 8, 10 };
		double[] by = new double[] { 0, 5, 10, 5, 0 };

		XYTrace traceB = new XYTrace(bx, by, GraphType.area, "name");
		
		assertEquals(0.782608687877655, XYTrace.correlate(traceA, traceB), 1e-12);

		ax = new double[] { 1, 3, 5, 7, 9 };
		ay = new double[] { 0, 5, 10, 5, 0 };

		traceA = new XYTrace(ax, ay, GraphType.area, "name");


		bx = new double[] { 2, 4, 6, 8, 10 };
		by = new double[] { 0, 3, 0, 5, 0 };

		traceB = new XYTrace(bx, by, GraphType.area, "name");
		
		assertEquals(0.11101597547531128, XYTrace.correlate(traceA, traceB), 1e-12);

		ax = new double[] { 1, 3, 5, 7, 9 };
		ay = new double[] { 0, 5, 10, 5, 0 };

		traceA = new XYTrace(ax, ay, GraphType.area, "name");


		bx = new double[] { 2, 4, 6, 8, 10 };
		by = new double[] { 2.5, 7.5, 7.5, 2.5, 0 };

		traceB = new XYTrace(bx, by, GraphType.area, "name");

		assertEquals(0.9548621773719788, XYTrace.correlate(traceA, traceB), 1e-12);

		ax = new double[] { 1, 3, 5, 7, 9 };
		ay = new double[] { 0, 5, 10, 5, 0 };

		traceA = new XYTrace(ax, ay, GraphType.area, "name");


		bx = new double[] { 1, 3, 5, 7, 9 };
		by = new double[] { 0, 3, 0, 5, 0 };

		traceB = new XYTrace(bx, by, GraphType.area, "name");
		
		assertEquals(0.20766964554786682, XYTrace.correlate(traceA, traceB), 1e-12);

		ax = new double[] { 1, 3, 5, 7, 9 };
		ay = new double[] { 0, 5, 10, 5, 0 };

		traceA = new XYTrace(ax, ay, GraphType.area, "name");


		bx = new double[] { 1, 3, 5, 7, 9 };
		by = new double[] { 0, 5, 10, 5, 0 };

		traceB = new XYTrace(bx, by, GraphType.area, "name");
		
		assertEquals(1, XYTrace.correlate(traceA, traceB), 1e-12);
	}
	
	public void testAlign() {

		double[] ax = new double[] { 1, 3, 5, 7, 9 };
		double[] ay = new double[] { 0, 5, 10, 5, 0 };

		XYTrace traceA = new XYTrace(ax, ay, GraphType.area, "name");


		double[] bx = new double[] { 2, 4, 6, 8, 10 };
		double[] by = new double[] { 0, 3, 0, 5, 0 };

		XYTrace traceB = new XYTrace(bx, by, GraphType.area, "name");
		
		ArrayList<XYPoint>[] values=XYTrace.alignXYPoints(traceA.getPoints(), traceB.getPoints());

        ArrayList<XYPoint> alignedA = values[0];
        ArrayList<XYPoint> alignedB = values[1];
        
        assertEquals("alignedA should have 10 points", 10, alignedA.size());
        assertEquals("alignedB should have 10 points", 10, alignedB.size());

        assertEquals("x=2 in alignedA", 2.0, alignedA.get(1).getX(), 1e-12);
        assertEquals("y(0) for A should be 0", 0, alignedA.get(0).getY(), 1e-12);
        assertEquals("y(1) for A should be 2.5", 2.5, alignedA.get(1).getY(), 1e-12);
        assertEquals("y(2) for A should be 5", 5, alignedA.get(2).getY(), 1e-12);

        assertEquals(0.0, alignedB.get(0).getY(), 1e-12);
        assertEquals(0.0, alignedB.get(1).getY(), 1e-12);
        assertEquals(1.5, alignedB.get(2).getY(), 1e-12);

		System.out.println(new XYTrace(values[0], GraphType.area, "name"));

		System.out.println(new XYTrace(values[1], GraphType.area, "name"));
	}
	
	/**
     * Test the behavior when both input lists are empty.
     */
    public void testBothEmpty() {
        ArrayList<XYPoint> listA = new ArrayList<>();
        ArrayList<XYPoint> listB = new ArrayList<>();
        
        ArrayList<XYPoint>[] result = XYTrace.alignXYPoints(listA, listB);
        ArrayList<XYPoint> alignedA = result[0];
        ArrayList<XYPoint> alignedB = result[1];
        
        assertEquals("Aligned A should be empty", 0, alignedA.size());
        assertEquals("Aligned B should be empty", 0, alignedB.size());
    }

    /**
     * Test the behavior when one list is empty and the other is not.
     */
    public void testOneEmpty() {
        ArrayList<XYPoint> listA = new ArrayList<>();
        listA.add(new XYPoint(1.0, 10.0));
        listA.add(new XYPoint(2.0, 20.0));
        
        ArrayList<XYPoint> listB = new ArrayList<>(); // empty

        ArrayList<XYPoint>[] result = XYTrace.alignXYPoints(listA, listB);
        ArrayList<XYPoint> alignedA = result[0];
        ArrayList<XYPoint> alignedB = result[1];

        // We expect every X from listA to appear in both, with B interpolated from null
        assertEquals("Aligned A should have same size as listA", 2, alignedA.size());
        assertEquals("Aligned B should also have 2 points", 2, alignedB.size());

        // Check that the X-values match and Y-values in B are from the interpolation rule
        assertEquals(1.0, alignedA.get(0).getX(), 1e-12);
        assertEquals(10.0, alignedA.get(0).getY(), 1e-12);
        
        assertEquals("Check Y from extrapolation rules", 0.0, alignedB.get(0).getY(), 1e-12);
    }

    /**
     * Test when both lists have exactly the same single point.
     */
    public void testSameSinglePoint() {
        ArrayList<XYPoint> listA = new ArrayList<>();
        listA.add(new XYPoint(0.0, 5.0));
        
        ArrayList<XYPoint> listB = new ArrayList<>();
        listB.add(new XYPoint(0.0, 5.0));

        ArrayList<XYPoint>[] result = XYTrace.alignXYPoints(listA, listB);
        ArrayList<XYPoint> alignedA = result[0];
        ArrayList<XYPoint> alignedB = result[1];
        
        assertEquals(1, alignedA.size());
        assertEquals(1, alignedB.size());

        assertEquals(0.0, alignedA.get(0).getX(), 1e-12);
        assertEquals(5.0, alignedA.get(0).getY(), 1e-12);
        assertEquals(0.0, alignedB.get(0).getX(), 1e-12);
        assertEquals(5.0, alignedB.get(0).getY(), 1e-12);
    }

    /**
     * Test partial overlap and interpolation.
     * For example: A has points at x=0 and x=2, B has points at x=1 and x=2.
     */
    public void testPartialOverlap() {
        // A: (0, 0) -> (2, 4)
        ArrayList<XYPoint> listA = new ArrayList<>();
        listA.add(new XYPoint(0.0, 0.0));
        listA.add(new XYPoint(2.0, 4.0));

        // B: (1, 10) -> (2, 8)
        ArrayList<XYPoint> listB = new ArrayList<>();
        listB.add(new XYPoint(1.0, 10.0));
        listB.add(new XYPoint(2.0, 8.0));

        ArrayList<XYPoint>[] result = XYTrace.alignXYPoints(listA, listB);
        ArrayList<XYPoint> alignedA = result[0];
        ArrayList<XYPoint> alignedB = result[1];

        // The code merges x=0,1,2. 
        // At x=0, B is interpolated from "prevB"=null or next B => see your code's rule. 
        // Then x=1 => direct from B, interpolate A
        // Then x=2 => direct from both.

        // We'll just check the final aligned lists are size=3
        assertEquals(3, alignedA.size());
        assertEquals(3, alignedB.size());

        // Check x-values in ascending order
        assertEquals(0.0, alignedA.get(0).getX(), 1e-12);
        assertEquals(1.0, alignedA.get(1).getX(), 1e-12);
        assertEquals(2.0, alignedA.get(2).getX(), 1e-12);

        assertEquals(0.0, alignedB.get(0).getX(), 1e-12);
        assertEquals(1.0, alignedB.get(1).getX(), 1e-12);
        assertEquals(2.0, alignedB.get(2).getX(), 1e-12);

        // You can also check the Y-values if you want to confirm the interpolation logic.
        // For example, at x=1 in listA, interpolation between (0,0) and (2,4) => y=2
        assertEquals("A's Y at x=1 should be halfway between 0 and 4", 2.0, alignedA.get(1).getY(), 1e-12);
        // B is direct at x=1 => y=10
        assertEquals(10.0, alignedB.get(1).getY(), 1e-12);
    }

    /**
     * Put everything together with a bigger test that checks multi-step interpolation.
     */
    public void testMultiStepInterpolation() {
        // A: 0->(1), 2->(5), 4->(9)
        ArrayList<XYPoint> listA = new ArrayList<>();
        listA.add(new XYPoint(0.0, 1.0));
        listA.add(new XYPoint(2.0, 5.0));
        listA.add(new XYPoint(4.0, 9.0));

        // B: 1->(3), 3->(7), 4->(10)
        ArrayList<XYPoint> listB = new ArrayList<>();
        listB.add(new XYPoint(1.0, 3.0));
        listB.add(new XYPoint(3.0, 7.0));
        listB.add(new XYPoint(4.0, 10.0));

        ArrayList<XYPoint>[] result = XYTrace.alignXYPoints(listA, listB);
        ArrayList<XYPoint> alignedA = result[0];
        ArrayList<XYPoint> alignedB = result[1];

        // We expect x=0,1,2,3,4 in the final alignment:
        assertEquals(5, alignedA.size());
        assertEquals(5, alignedB.size());

        // Check the ascending X
        assertEquals(0.0, alignedA.get(0).getX(), 1e-12);
        assertEquals(1.0, alignedA.get(1).getX(), 1e-12);
        assertEquals(2.0, alignedA.get(2).getX(), 1e-12);
        assertEquals(3.0, alignedA.get(3).getX(), 1e-12);
        assertEquals(4.0, alignedA.get(4).getX(), 1e-12);

        // And B's X as well
        assertEquals(0.0, alignedB.get(0).getX(), 1e-12);
        assertEquals(1.0, alignedB.get(1).getX(), 1e-12);
        assertEquals(2.0, alignedB.get(2).getX(), 1e-12);
        assertEquals(3.0, alignedB.get(3).getX(), 1e-12);
        assertEquals(4.0, alignedB.get(4).getX(), 1e-12);

        // Spot check Y-values as needed...
    }
}
