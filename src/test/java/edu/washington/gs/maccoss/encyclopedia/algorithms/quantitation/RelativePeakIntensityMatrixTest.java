package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import junit.framework.TestCase;

public class RelativePeakIntensityMatrixTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false, true, false);
	private static final AminoAcidConstants aaConstants=PARAMETERS.getAAConstants();
	
	public void testELHPVLK() {
		double[] massArray=new double[] {130.0498694668129, 147.1128041505165, 243.1339334668129, 260.1968681505165, 359.2652821505165, 380.19284546681286, 456.31804615051647, 576.3140234668128, 593.3769581505164};float[] intensityArray=new float[] {773459.44f, 5570774.0f, 376921.66f, 3897037.5f, 605428.9f, 2419830.0f, 5553414.5f, 13243.6455f, 268873.78f};float[] correlationArray=new float[] {-0.37315476f, 0.9890311f, 0.79180455f, 0.99583536f, 0.99454385f, 0.9989461f, 0.9980313f, 0.70449966f, 0.8637193f};
		LibraryEntry XXX_2019_0304_RJ_12_2_ELHPVLK=new LibraryEntry("XXX_2019_0304_RJ_12_2.mzML", new HashSet<String>(), 1, 418.2554458086647, (byte)2, "ELHPVLK", 1, 2835.5513f, 0.0930104f, massArray, intensityArray, correlationArray, aaConstants);

		massArray=new double[] {130.0498694668129, 147.1128041505165, 243.1339334668129, 260.1968681505165, 359.2652821505165, 380.19284546681286, 456.31804615051647, 477.2456094668129, 576.3140234668128, 593.3769581505164, 689.3980874668129};intensityArray=new float[] {159851.83f, 7726639.0f, 582203.8f, 5575229.5f, 715657.75f, 3258831.5f, 7283809.5f, 143028.08f, 187806.05f, 625733.7f, 16991.264f};correlationArray=new float[] {-0.18069786f, 0.97763574f, 0.91893303f, 0.9945558f, 0.9075433f, 0.9971901f, 0.9974563f, 0.8843416f, 0.8473087f, 0.9847535f, 0.5808859f};
		LibraryEntry XXX_2019_0304_RJ_13_3_ELHPVLK=new LibraryEntry("XXX_2019_0304_RJ_13_3.mzML", new HashSet<String>(), 1, 418.2554458086647, (byte)2, "ELHPVLK", 1, 2852.5515f, 0.133916f, massArray, intensityArray, correlationArray, aaConstants);

		massArray=new double[] {130.0498694668129, 147.1128041505165, 243.1339334668129, 260.1968681505165, 359.2652821505165, 380.19284546681286, 456.31804615051647, 477.2456094668129, 576.3140234668128, 593.3769581505164, 689.3980874668129};intensityArray=new float[] {163603.72f, 3790984.0f, 222426.69f, 2890778.2f, 462041.16f, 1987162.5f, 4051529.0f, 38186.75f, 38108.113f, 255911.22f, 14044.433f};correlationArray=new float[] {0.89856327f, 0.99607205f, 0.9355579f, 0.9970543f, 0.9729357f, 0.99407506f, 0.99596477f, 0.7329936f, 0.85431474f, 0.9419776f, 0.92073727f};
		LibraryEntry XXX_2019_0304_RJ_16_4_ELHPVLK=new LibraryEntry("XXX_2019_0304_RJ_16_4.mzML", new HashSet<String>(), 1, 418.2554458086647, (byte)2, "ELHPVLK", 1, 2863.515f, 0.124273f, massArray, intensityArray, correlationArray, aaConstants);

		massArray=new double[] {130.0498694668129, 147.1128041505165, 243.1339334668129, 260.1968681505165, 359.2652821505165, 380.19284546681286, 456.31804615051647, 576.3140234668128, 593.3769581505164};intensityArray=new float[] {14467.836f, 880617.75f, 95044.91f, 515502.16f, 94848.92f, 318728.44f, 642752.25f, 49509.527f, 123847.586f};correlationArray=new float[] {0.41985622f, 0.9981061f, 0.6985816f, 0.99814504f, 0.51362395f, 0.9957508f, 0.9683793f, 1.4E-45f, 0.9110258f};
		LibraryEntry XXX_2019_0304_RJ_23_9_ELHPVLK=new LibraryEntry("XXX_2019_0304_RJ_23_9.mzML", new HashSet<String>(), 1, 418.2554458086647, (byte)2, "ELHPVLK", 1, 2850.8184f, 0.126865f, massArray, intensityArray, correlationArray, aaConstants);

		massArray=new double[] {130.0498694668129, 147.1128041505165, 243.1339334668129, 260.1968681505165, 359.2652821505165, 380.19284546681286, 456.31804615051647, 477.2456094668129, 576.3140234668128, 593.3769581505164, 689.3980874668129};intensityArray=new float[] {104980.695f, 4103354.0f, 343504.62f, 3175719.0f, 471142.25f, 2018141.4f, 4500639.5f, 25821.865f, 71118.766f, 220849.12f, 16920.418f};correlationArray=new float[] {-0.10478447f, 0.9896324f, 0.9108415f, 0.9969171f, 0.9948849f, 0.996705f, 0.9969578f, 0.97556734f, 0.9302268f, 0.93036294f, 0.97556734f};
		LibraryEntry XXX_2019_0304_RJ_24_12_ELHPVLK=new LibraryEntry("XXX_2019_0304_RJ_24_12.mzML", new HashSet<String>(), 1, 418.2554458086647, (byte)2, "ELHPVLK", 1, 2797.564f, 0.0887852f, massArray, intensityArray, correlationArray, aaConstants);

		massArray=new double[] {130.0498694668129, 147.1128041505165, 243.1339334668129, 260.1968681505165, 359.2652821505165, 380.19284546681286, 456.31804615051647, 477.2456094668129, 576.3140234668128, 593.3769581505164, 689.3980874668129};intensityArray=new float[] {65053.9f, 1.0289069E7f, 833442.2f, 8275591.0f, 1476940.1f, 5438929.0f, 1.118502E7f, 272807.78f, 504273.78f, 915608.06f, 190840.88f};correlationArray=new float[] {-0.2768532f, 0.9981788f, 0.9839694f, 0.99954575f, 0.99923795f, 0.9991896f, 0.9995219f, 0.9710275f, 0.97830456f, 0.99317473f, 0.99216247f};
		LibraryEntry XXX_2019_0304_RJ_26_14_ELHPVLK=new LibraryEntry("XXX_2019_0304_RJ_26_14.mzML", new HashSet<String>(), 1, 418.2554458086647, (byte)2, "ELHPVLK", 1, 2818.936f, 0.0736802f, massArray, intensityArray, correlationArray, aaConstants);

		massArray=new double[] {130.0498694668129, 147.1128041505165, 243.1339334668129, 260.1968681505165, 359.2652821505165, 380.19284546681286, 456.31804615051647, 593.3769581505164};intensityArray=new float[] {128095.63f, 3343094.5f, 278725.7f, 2435406.0f, 391872.34f, 1397072.1f, 3151301.0f, 217233.53f};correlationArray=new float[] {0.98320955f, 0.99541736f, 0.2451348f, 0.99387455f, 0.99803084f, 0.9990833f, 0.9995974f, 0.9474021f};
		LibraryEntry XXX_2019_0304_RJ_34_20_ELHPVLK=new LibraryEntry("XXX_2019_0304_RJ_34_20.mzML", new HashSet<String>(), 1, 418.2554458086647, (byte)2, "ELHPVLK", 1, 2783.1426f, 0.113011f, massArray, intensityArray, correlationArray, aaConstants);

		RelativePeakIntensityMatrix matrix=new RelativePeakIntensityMatrix(PARAMETERS.getFragmentTolerance(), Optional.empty());
		matrix.addPeak("XXX_2019_0304_RJ_12_2", XXX_2019_0304_RJ_12_2_ELHPVLK);
		matrix.addPeak("XXX_2019_0304_RJ_13_3", XXX_2019_0304_RJ_13_3_ELHPVLK);
		matrix.addPeak("XXX_2019_0304_RJ_16_4", XXX_2019_0304_RJ_16_4_ELHPVLK);
		matrix.addPeak("XXX_2019_0304_RJ_23_9", XXX_2019_0304_RJ_23_9_ELHPVLK);
		matrix.addPeak("XXX_2019_0304_RJ_24_12", XXX_2019_0304_RJ_24_12_ELHPVLK);
		matrix.addPeak("XXX_2019_0304_RJ_26_14", XXX_2019_0304_RJ_26_14_ELHPVLK);
		matrix.addPeak("XXX_2019_0304_RJ_34_20", XXX_2019_0304_RJ_34_20_ELHPVLK);
		
		double[] targets=new double[] {260.1968681505165, 147.1128041505165, 456.31804615051647, 380.19284546681286,
				359.2652821505165, 243.1339334668129, 593.3769581505164, 576.3140234668128, 477.2456094668129,
				130.0498694668129, 689.3980874668129};
		
		for (int n = 1; n <= targets.length; n++) {
			double[] masses=matrix.pickNBestPeaks(n, 0);
			assertEquals(Math.min(7, n), masses.length);
			for (int i = 0; i < Math.min(7, n); i++) {
				assertEquals(targets[i], masses[i], 0.0001);
			}
		}
	}
}
