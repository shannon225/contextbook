package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.awt.Dimension;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import junit.framework.TestCase;

public class LibraryEntryTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"));
	public void testReverse() {
		double[] massArray = new double[] { 98.06063, 175.11955, 227.10323,
				304.16214, 324.15599, 333, 419.18908, 444, 505.20367,
				532.27314, 555, 618.28773, 650.407259, 666, 713.32082,
				733.31467, 777, 779.449849, 810.37359, 862.35727, 876.502609,
				888, 939.41618, 1018.45838, 1036.46894 };
		float[] intensityArray = new float[] { 1f, 2f, 3f, 4f, 5f, 6f,
				7f, 8f, 9f, 10f, 11f, 12f, 13f, 14f, 15f, 16f, 17f, 18f, 19f,
				20f, 21f, 22f, 23f, 24f, 25f };
		
		LibraryEntry entry=new LibraryEntry("", new HashSet<String>(), 518.73841, (byte)2, "PEPT[+80]IDER", 1, 0.0f, 0.0f, massArray, intensityArray);
		LibraryEntry reverse=entry.getDecoy(PARAMETERS);
		assertEquals("EDIT[+79.966331]PEPR", reverse.getPeptideModSeq());
		
		double[] reverseMasses=reverse.getMassArray();
		double[] expectedReversedMasses=new double[] { 130.05045900000002, 175.11955, 245.07740900000002, 272.172311, 333.0, 358.161469, 401.21490099999994, 444.0, 498.267661, 505.20367,
				555.0, 618.28773, 650.407259, 666.0, 713.32082, 733.31467, 777.0, 779.449849, 810.37359, 862.3572700000001, 876.502609, 888.0, 939.41618, 1018.4583800000003,
				1036.4689399999997 };
		for (int i = 0; i < reverseMasses.length; i++) {
			assertEquals(expectedReversedMasses[i], reverseMasses[i], 0.1);
		}
	}
	public static void main(String[] args) {
		HashMap<String, String> defaults=SearchParameterParser.getDefaultParameters();
		defaults.put("-ptol", "17");
		defaults.put("-ftol", "17");
		defaults.put("-lftol", "17");
		SearchParameters parameters=SearchParameterParser.parseParameters(defaults);
		
		double[] massArray=new double[] { 161.32635498046875, 170.09095764160156, 171.07424926757812, 175.1168975830078, 175.12811279296875, 179.0912628173828, 181.09426879882812, 185.56866455078125,
				188.1009063720703, 189.1050567626953, 194.31455993652344, 199.10459899902344, 207.0857391357422, 216.09510803222656, 225.09564208984375, 236.49539184570312, 242.58338928222656,
				244.3992462158203, 248.4324493408203, 258.7773742675781, 269.1279602050781, 278.11016845703125, 295.1349792480469, 304.1467590332031, 313.1484375, 335.1434326171875,
				360.16326904296875, 372.15997314453125, 401.209716796875, 402.2124328613281, 408.6524353027344, 417.9472961425781, 450.20111083984375, 470.2284851074219, 501.7383117675781,
				502.5447998046875, 504.5861511230469, 508.2248229980469, 531.5833129882812, 536.9232177734375, 537.25390625, 537.5802612304688, 537.9183959960938, 543.2594604492188, 563.9029541015625,
				564.9290771484375, 569.5771484375, 569.9093017578125, 570.2459716796875, 570.5816650390625, 580.27490234375, 580.7698364257812, 585.7479248046875, 591.9186401367188, 592.2725830078125,
				592.7791137695312, 597.2657470703125, 602.9380493164062, 603.2718505859375, 603.6080322265625, 604.769287109375, 607.287841796875, 608.2918701171875, 609.2989501953125,
				635.2683715820312, 635.6048583984375, 636.2701416015625, 640.7933959960938, 641.2754516601562, 648.8037109375, 649.3048095703125, 649.8116455078125, 674.2930297851562,
				703.828369140625, 704.31298828125, 704.823974609375, 717.3251953125, 718.3211059570312, 735.3412475585938, 736.3433837890625, 753.3282470703125, 756.3518676757812, 756.8523559570312,
				761.3369140625, 770.8458251953125, 795.869873046875, 796.3707275390625, 796.8701171875, 797.3671875, 804.8756713867188, 805.3729858398438, 805.8749389648438, 806.3773803710938,
				806.86328125, 813.880615234375, 814.379150390625, 844.8610229492188, 845.3523559570312, 845.8623657226562, 848.4252319335938, 853.8621215820312, 854.3611450195312, 854.8626708984375,
				855.3735961914062, 855.89501953125, 874.4122924804688, 935.4593505859375, 936.4617919921875, 1002.4571533203125, 1034.4339599609375, 1072.511474609375, 1073.519775390625,
				1090.5216064453125, 1159.5457763671875, 1208.546875 };
		float[] intensityArray=new float[] { 789.0159f, 1386.6464f, 15342.616f, 9691.694f, 731.12445f, 1063.6119f, 1085.3112f, 1121.8958f, 22498.57f, 2170.6208f, 841.44977f, 991.2484f, 1041.3636f,
				11757.499f, 4665.62f, 958.1826f, 953.03204f, 885.48535f, 1082.2528f, 1037.7479f, 1041.7932f, 1003.23114f, 2584.1646f, 1223.9075f, 2429.5198f, 1669.2076f, 1105.715f, 1264.8138f,
				27486.086f, 5866.4863f, 1360.7203f, 1256.7516f, 2053.1938f, 7953.5073f, 1525.7118f, 1914.5751f, 1217.0173f, 1885.7172f, 1708.9697f, 11138.359f, 12531.536f, 6699.1206f, 1757.6083f,
				1407.471f, 1838.0353f, 1537.4362f, 15354.347f, 13217.588f, 1987.985f, 1521.9126f, 4462.602f, 1763.949f, 1694.7002f, 2851.2668f, 1820.5972f, 1965.2032f, 7103.8667f, 4117.0996f,
				4849.963f, 5098.8604f, 2546.4956f, 15416.269f, 6805.201f, 2060.623f, 2078.2737f, 4805.7695f, 2077.287f, 2386.0088f, 1427.156f, 7512.233f, 4608.3545f, 2036.2054f, 4946.5f, 1852.7395f,
				4057.4482f, 2282.1423f, 1639.9491f, 1750.5055f, 8599.824f, 1871.3995f, 4476.562f, 5776.7046f, 5573.211f, 5098.2173f, 2677.7788f, 7161.399f, 11230.046f, 8512.083f, 1711.4738f,
				42050.598f, 36705.426f, 21339.602f, 6402.2563f, 1767.8688f, 2210.962f, 6152.833f, 1921.9388f, 2189.6343f, 1838.1809f, 5342.611f, 6875.591f, 14617.56f, 2743.2551f, 1799.31f, 1628.9028f,
				1726.5023f, 7531.061f, 2296.9006f, 2204.0225f, 2103.5925f, 4470.2554f, 2748.8965f, 1818.1165f, 2741.7195f, 2618.4417f };
		float[] correlationArray=new float[] { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
				1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
				1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
				1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
		LibraryEntry entry=new LibraryEntry("VillenJ_Exactive_HumanPhosphoproteome.blib", new HashSet<String>(), 1, 641.280032, (byte) 3, "NTPS[+79.966331]QHSHSIQHSPER", 2, 33.986168f, 2.04E-4f, massArray, intensityArray,
				correlationArray);
		
		AnnotatedLibraryEntry annotated=new AnnotatedLibraryEntry(entry, parameters);
		for (int ion=0; ion<annotated.getIonAnnotations().length; ion++) {
			FragmentIon fragmentIon=annotated.getIonAnnotations()[ion];
			if (fragmentIon!=null) {
				System.out.println(ion+") "+fragmentIon.toCanonicalIonTypeString());
			}
		}

		Charter.launchComponent(Charter.getChart(annotated), annotated.getName(), new Dimension(1000, 500));
	}
}
