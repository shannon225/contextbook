package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.Random;

import junit.framework.TestCase;

public class SGFilterTest extends TestCase {

	public void testSavitzkyGolaySmooth() {
		float[] data=new float[] {78122.11f, 98378.266f, 142467.98f, 160690.98f, 222905.19f, 232847.81f, 307916.06f, 287985.97f, 354645.78f, 292500.47f, 363600.66f, 315389.8f, 347600.38f, 301699.44f,
				294894.66f, 221306.75f, 186132.05f, 170257.3f, 146349.67f};
		float[] expected=new float[] {78122.11f, 98378.266f, 142467.98f, 160690.98f, 209040.94f, 253904.34f, 279367.1f, 309400.34f, 321014.0f, 335611.34f, 334470.47f, 340195.16f, 326188.44f,
				309110.4f, 270765.3f, 221306.75f, 186132.05f, 170257.3f, 146349.67f};

		float[] smooth=SkylineSGFilter.savitzkyGolaySmooth(data);
		for (int i=0; i<smooth.length; i++) {
			assertEquals(expected[i], smooth[i], 1);
		}

		data=new float[] {3496.0784f, 3709.9358f, 9963.435f, 7364.1235f, 14094.167f, 10891.73f, 20391.303f, 6583.6714f, 20247.096f, 19192.396f, 17108.42f, 12707.958f, 21243.793f, 9334.771f,
				15843.864f, 3872.4568f, 11016.68f, 4989.8164f, 7633.1284f};
		expected=new float[] {3496.0784f, 3709.9358f, 9963.435f, 7364.1235f, 11457.621f, 12947.143f, 14239.61f, 16326.898f, 15745.773f, 17576.006f, 16502.514f, 18264.258f, 14697.832f, 13359.507f,
				11404.038f, 3872.4568f, 11016.68f, 4989.8164f, 7633.1284f};

		smooth=SkylineSGFilter.savitzkyGolaySmooth(data);
		for (int i=0; i<smooth.length; i++) {
			assertEquals(expected[i], smooth[i], 1);
		}

		data=new float[] {0f, 0f, 0f, 2272.649f, 2853.767f, 2879.7114f, 2626.8823f, 0f, 3799.889f, 7058.896f, 9528.616f, 0f, 2981.3376f, 1523.2601f, 6134.175f, 955.7171f, 0f, 0f, 2771.0017f};
		expected=new float[] {0.0f, 0.0f, 0.0f, 2272.649f, 2031.3876f, 1988.9738f, 2166.7698f, 3724.152f, 4316.957f, 4858.2983f, 4524.589f, 4888.5215f, 3902.211f, 2617.0857f, 1531.434f, 955.7171f,
				0.0f, 0.0f, 2771.0017f};
		smooth=SkylineSGFilter.savitzkyGolaySmooth(data);
		for (int i=0; i<smooth.length; i++) {
			assertEquals(expected[i], smooth[i], 1);
		}
	}


	public void testIterativePaddedSavitzkyGolaySmooth() {
		Random rd = new Random();
		for (int n = 0; n < 10000; n++) {
			float[] arr = new float[1000];
			for (int i = 0; i < arr.length; i++) {
				arr[i] = rd.nextFloat();
			}

			float[] smooth=SkylineSGFilter.paddedSavitzkyGolaySmooth(arr);
			assertTrue(smooth.length==arr.length);
		}
	}

	public void testPaddedSavitzkyGolaySmooth() {
		float[] data=new float[] {78122.11f, 98378.266f, 142467.98f, 160690.98f, 222905.19f, 232847.81f, 307916.06f, 287985.97f, 354645.78f, 292500.47f, 363600.66f, 315389.8f, 347600.38f, 301699.44f,
				294894.66f, 221306.75f, 186132.05f, 170257.3f, 146349.67f};
		float[] expected=new float[] {56478.56f, 96164.49f, 133891.97f, 179591.2f, 209040.94f, 253904.34f, 279367.1f, 309400.34f, 321014.0f, 335611.34f, 334470.47f, 340195.16f, 326188.44f, 309110.4f,
				270765.3f, 249917.61f, 200255.25f, 149017.28f, 95208.57f};

		float[] smooth=SkylineSGFilter.paddedSavitzkyGolaySmooth(data);
		for (int i=0; i<smooth.length; i++) {
			assertEquals(expected[i], smooth[i], 1);
		}

		data=new float[] {3496.0784f, 3709.9358f, 9963.435f, 7364.1235f, 14094.167f, 10891.73f, 20391.303f, 6583.6714f, 20247.096f, 19192.396f, 17108.42f, 12707.958f, 21243.793f, 9334.771f,
				15843.864f, 3872.4568f, 11016.68f, 4989.8164f, 7633.1284f};
		expected=new float[] {2607.3567f, 5201.267f, 6909.646f, 10819.156f, 11457.621f, 12947.143f, 14239.61f, 16326.898f, 15745.773f, 17576.006f, 16502.514f, 18264.258f, 14697.832f, 13359.507f,
				11404.038f, 10281.441f, 7483.6216f, 6399.561f, 3770.3394f};

		smooth=SkylineSGFilter.paddedSavitzkyGolaySmooth(data);
		for (int i=0; i<smooth.length; i++) {
			assertEquals(expected[i], smooth[i], 1);
		}
		System.out.println();

		data=new float[] {0f, 0f, 0f, 2272.649f, 2853.767f, 2879.7114f, 2626.8823f, 0f, 3799.889f, 7058.896f, 9528.616f, 0f, 2981.3376f, 1523.2601f, 6134.175f, 955.7171f, 0f, 0f, 2771.0017f};
		expected=new float[] {0.0f, 294.8576f, 948.7938f, 1892.9645f, 2031.3876f, 1988.9738f, 2166.7698f, 3724.152f, 4316.957f, 4858.2983f, 4524.589f, 4888.5215f, 3902.211f, 2617.0857f, 1531.434f,
				2283.8643f, 1548.174f, 1042.4114f, 208.01495f};
		smooth=SkylineSGFilter.paddedSavitzkyGolaySmooth(data);
		for (int i=0; i<smooth.length; i++) {
			assertEquals(expected[i], smooth[i], 1);
		}
	}
}
