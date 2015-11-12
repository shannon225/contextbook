package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.Arrays;

import junit.framework.TestCase;

public class FragmentationModelTest extends TestCase {
	public void testPrimaryIons() {
		String sequence="PEPT[+80]IDER";
		FragmentationModel model=new FragmentationModel(sequence);
		double[] ions=model.getPrimaryIons(FragmentationType.CID);
		double[] expected=new double[] {98.06063, 175.11955, 227.10323, 304.16215, 324.15603, 407.226834, 419.18915, 505.20373, 520.310934, 532.27325, 615.344054, 618.28783, 635.337934,
				712.3968540000001, 713.32095, 733.31483, 764.380534, 810.3737500000001, 841.4394540000001, 862.35743, 920.481634, 938.4922540000001, 939.4163500000001, 1018.45853, 1036.4691500000001};
		for (int i=0; i<ions.length; i++) {
			assertEquals(expected[i], ions[i], 0.0001);
		}
		
		ions=model.getPrimaryIons(FragmentationType.ETD);
		expected=new double[] {115.08717911000001, 158.09300089, 159.10082593, 244.12977911, 287.13560089, 288.14342593, 341.18257911, 402.16260088999996, 403.17042592999996, 424.25338311,
				515.2467008899999, 516.25452593, 522.23027911, 537.33748311, 598.31750489, 599.3253299300001, 635.31437911, 652.36448311, 695.37030489, 696.29440089, 696.3781299300001,
				697.3022259300001, 750.34137911, 781.40708311, 793.3472008900001, 794.3550259300001, 824.41290489, 825.4207299300001, 879.38397911, 921.4657048900001, 922.3898008900001,
				922.4735299300002, 923.3976259300001, 937.50818311, 1019.4426008900001, 1020.4504259300002, 1035.48507911};
		for (int i=0; i<ions.length; i++) {
			assertEquals(expected[i], ions[i], 0.0001);
		}
	}
	
	public void testFragmentation() {
		String sequence="PEPTIDER";
		double[] expectedB = new double[] { 98.06063, 227.10323, 324.15599, 425.20367, 538.28773, 653.31467, 782.35727, 938.45838 };
		double[] expectedY = new double[] { 175.11955, 304.16214, 419.18908, 532.27314, 633.32082, 730.37359, 859.41618, 956.46894 };
		
		FragmentationModel model=new FragmentationModel(sequence);
		double[] bs=model.getBIons();
		for (int i = 0; i < bs.length; i++) {
			assertEquals(expectedB[i], bs[i], 0.001);
		}
		double[] ys=model.getYIons();
		for (int i = 0; i < ys.length; i++) {
			assertEquals(expectedY[i], ys[i], 0.001);
		}
	}
	public void testNeutralLossFragmentation() {
		String sequence="PEPT[+80]IDER";
		double[] expectedB = new double[] { 98.06063, 227.10323, 324.15599, 425.20367+80, 538.28773+80, 653.31467+80, 782.35727+80, 938.45838+80, 425.20367+80-97.976896, 538.28773+80-97.976896, 653.31467+80-97.976896, 782.35727+80-97.976896, 938.45838+80-97.976896  };
		Arrays.sort(expectedB);
		double[] expectedY = new double[] { 175.11955, 304.16214, 419.18908, 532.27314, 633.32082+80, 730.37359+80, 859.41618+80, 956.46894+80, 633.32082+80-97.976896, 730.37359+80-97.976896, 859.41618+80-97.976896, 956.46894+80-97.976896  };
		Arrays.sort(expectedY);
		
		FragmentationModel model=new FragmentationModel(sequence);
		double[] bs=model.getBIons();
		for (int i = 0; i < bs.length; i++) {
			assertEquals(expectedB[i], bs[i], 0.001);
		}
		double[] ys=model.getYIons();
		for (int i = 0; i < ys.length; i++) {
			assertEquals(expectedY[i], ys[i], 0.001);
		}
	}
	
	public void testGetMasses() {
		String sequence="PEPTIDER";
		double[] expected = new double[] { 97.0528, 129.0426, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011 };
		double[] masses=FragmentationModel.getMasses(sequence).x;
		for (int i = 0; i < masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.0001);
		}
		
		sequence="PEPT[+80]IDER";
		expected = new double[] { 97.0528, 129.0426, 97.0528, 101.0477+80.0, 113.0841, 115.027, 129.0426, 156.1011 };
		masses=FragmentationModel.getMasses(sequence).x;
		for (int i = 0; i < masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.0001);
		}
		
		sequence="PE[-17]PTIDER";
		expected = new double[] { 97.0528, 129.0426-17.0, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011 };
		masses=FragmentationModel.getMasses(sequence).x;
		for (int i = 0; i < masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.0001);
		}
		
		sequence="[-17]PEPTIDER";
		expected = new double[] { 97.0528-17.0, 129.0426, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011 };
		masses=FragmentationModel.getMasses(sequence).x;
		for (int i = 0; i < masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.0001);
		}
		
		sequence="[+42]PEPTIDER";
		expected = new double[] { 97.0528+42.0, 129.0426, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011 };
		masses=FragmentationModel.getMasses(sequence).x;
		for (int i = 0; i < masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.0001);
		}
		
		sequence="PEPTIDER[+14]";
		expected = new double[] { 97.0528, 129.0426, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011+14.0 };
		masses=FragmentationModel.getMasses(sequence).x;
		for (int i = 0; i < masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.0001);
		}
	}
}
