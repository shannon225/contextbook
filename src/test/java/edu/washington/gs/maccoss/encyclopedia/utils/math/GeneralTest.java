package edu.washington.gs.maccoss.encyclopedia.utils.math;

import edu.washington.gs.maccoss.encyclopedia.datastructures.IntRange;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import junit.framework.TestCase;

public class GeneralTest extends TestCase {
	public static void main(String[] args) {
		float[] trace1=new float[] {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,112.29789f,423.9019f,1224.1725f,2097.833f,2401.0515f,2442.4944f,1957.5907f,1112.6595f,663.86957f,70.27865f,90.12375f,452.2014f,778.6451f,1336.3347f,1713.0786f,1695.5273f,1582.425f,1025.2372f};
		float[] trace2=new float[] {220.44585f,240.8575f,220.44585f,159.21089f,57.152626f,0.0f,0.0f,0.0f,0.0f,39.170887f,109.11889f,151.08769f,165.0773f,151.08769f,0.0f,0.0f,389.32227f,1019.15393f,1640.3691f,1752.7599f,2154.8984f,2047.731f,1431.2571f,746.2825f,748.64935f};
		float[] trace3=new float[] {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,23.815647f,66.34359f,91.86036f,33.17489f,136.65439f,191.12698f,125.09421f,200.71704f,305.5598f,308.63654f,245.67061f,116.66209f,132.78282f,47.665627f,0.0f,0.0f};
		float[] trace4=new float[] {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,255.44011f,633.43134f,1312.1504f,1667.1926f,1698.5582f,1406.2466f,790.25836f,189.12079f,16.356602f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
		float[] trace5=new float[] {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,58.54496f,42.353264f,306.30713f,350.11023f,616.84064f,643.0764f,735.50006f,631.24084f,606.2967f,338.6701f,295.63544f,34.493084f,55.7604f,0.0f,0.0f,0.0f,0.0f,0.0f};
		float[] rts=new float[trace1.length];
		for (int i=0; i<rts.length; i++) {
			rts[i]=i;
		}
		IntRange range=new IntRange(8, 18);
		
		if (false) {
			trace1=General.normalize(trace1, range);
			trace2=General.normalize(trace2, range);
			trace3=General.normalize(trace3, range);
			trace4=General.normalize(trace4, range);
			trace5=General.normalize(trace5, range);
		} else {
			trace1=General.normalizeAndBackgroundSubtract(trace1, range);
			trace2=General.normalizeAndBackgroundSubtract(trace2, range);
			trace3=General.normalizeAndBackgroundSubtract(trace3, range);
			trace4=General.normalizeAndBackgroundSubtract(trace4, range);
			trace5=General.normalizeAndBackgroundSubtract(trace5, range);
		}

		Charter.launchChart("Intens", "RT", false, getTrace(trace1, rts), getTrace(trace2, rts), getTrace(trace3, rts), getTrace(trace4, rts), getTrace(trace5, rts));
	}
	
	public void testNormalization() {
		float[] f=new float[10];
		for (int i = 0; i < 100; i++) {
			f[5]=100.0f*(float)Math.random();
			assertEquals(1.0f, General.normalizeToL2(f)[5]);
		}
	}
	
	public void testNumberOfOccurances() {
		String s="THISISATHESISTEST";
		assertEquals(2, General.numberOfOccurances(s, "TH"));
		assertEquals(3, General.numberOfOccurances(s, "IS"));
	}
	
	public void testConcatenate() {
		float[] a1=new float[] {1, 2, 3};
		float[] a2=new float[] {4, 5, 6};
		float[] a3=new float[] {7, 8, 9};

		float[] r=General.concatenate(a1, a2, a3);
		float[] expected=new float[] {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f};
		for (int i=0; i<r.length; i++) {
			assertEquals(expected[i], r[i]);
		}
	}
	
	public void testInsert() {
		String[] s1=new String[] {"a", "b", "c", "d", "e"};
		String[] s2=new String[] {"1", "2", "3"};
		String[] r=General.concatenate(s1, s2);
		String[] expected=new String[] {"a", "b", "c", "d", "e", "1", "2", "3"};
		for (int i = 0; i < r.length; i++) {
			assertEquals(r[i], expected[i]);
		}
		
		r=General.insert(s1, 0, s2);
		expected=new String[] {"1", "2", "3", "a", "b", "c", "d", "e"};
		for (int i = 0; i < r.length; i++) {
			assertEquals(r[i], expected[i]);
		}

		r=General.insert(s1, 2, s2);
		expected=new String[] {"a", "b", "1", "2", "3", "c", "d", "e"};
		for (int i = 0; i < r.length; i++) {
			assertEquals(r[i], expected[i]);
		}

		r=General.insert(s1, 5, s2);
		expected=new String[] {"a", "b", "c", "d", "e", "1", "2", "3"};
		for (int i = 0; i < r.length; i++) {
			assertEquals(r[i], expected[i]);
		}
	}
	
	public void testArrayIncrement() {
		int[] indices=new int[2];
		for (int i=0; i<1000; i++) {
			indices[0]++;
			indices[1]++;
		}
		assertEquals("1000,1000", General.toString(indices));
	}

	public static void main2(String[] args) {
	//public void testNormalize() {
		float[] b3=new float[] { 23338.361328125f, 16978.677734375f, 26238.6640625f, 28618.11328125f, 47211.97265625f, 60493.10546875f, 85625.6953125f, 154640.59375f, 163637.515625f, 113405.609375f,
				164475.375f, 202257.890625f, 100290.7734375f, 63675.58984375f, 31520.583984375f, 22526.6953125f, 0.0f, 0.0f, 0.0f, 6942.896484375f, 25359.82421875f, 26355.232421875f, 28414.279296875f,
				32256.48046875f, 28046.2421875f };
		
		float[] rts=new float[] { 30.342016220092773f, 30.380882263183594f, 30.42345428466797f, 30.462739944458008f, 30.503339767456055f, 30.543596267700195f, 30.583803176879883f, 30.622554779052734f,
				30.664770126342773f, 30.703386306762695f, 30.74576187133789f, 30.78369140625f, 30.825769424438477f, 30.865869522094727f, 30.906150817871094f, 30.945362091064453f, 30.98526954650879f,
				31.024911880493164f, 31.067047119140625f, 31.105043411254883f, 31.147865295410156f, 31.185548782348633f, 31.228477478027344f, 31.267902374267578f, 31.308513641357422f };
		
		float[] trace1=General.normalize(b3);
		float[] trace2=General.normalize(b3, new IntRange(3, b3.length-6));
		float[] trace3=General.normalizeAndBackgroundSubtract(b3, new IntRange(3, b3.length-6));
		float[] trace4=General.normalizeAndBackgroundSubtract(b3, new IntRange(0, b3.length-1));
		Charter.launchChart("Intens", "RT", false, getTrace(trace1, rts), getTrace(trace2, rts), getTrace(trace3, rts), getTrace(trace4, rts));
	}

	private static XYTraceInterface getTrace(float[] b3, float[] rts) {
		XYTraceInterface xyTrace=new XYTrace(General.toDoubleArray(rts), General.toDoubleArray(b3), GraphType.line, "B3");
		return xyTrace;
	}
}
