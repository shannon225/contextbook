package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.Arrays;

import junit.framework.TestCase;

public class FragmentationModelTest extends TestCase {
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
