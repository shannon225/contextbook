package edu.washington.gs.maccoss.encyclopedia.utils.math;

import edu.washington.gs.maccoss.encyclopedia.datastructures.IntRange;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import junit.framework.TestCase;

public class GeneralTest extends TestCase {
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

	public static void main(String[] args) {
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

	private static XYTrace getTrace(float[] b3, float[] rts) {
		XYTrace xyTrace=new XYTrace(General.toDoubleArray(rts), General.toDoubleArray(b3), GraphType.line, "B3");
		return xyTrace;
	}
}
