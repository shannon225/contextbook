package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;
import gnu.trove.list.array.TDoubleArrayList;

public class TransitionRefinerTest {
	public static void main(String[] args) {
		ArrayList<float[]> chromatograms=new ArrayList<float[]>();
		float[] y2=new float[] {37789.328125f, 37134.6796875f, 23441.00390625f, 41935.12890625f, 61996.5234375f, 46208.08984375f, 42318.9375f, 84970.9765625f, 105310.4375f, 183981.203125f,
				212268.0625f, 284312.4375f, 168349.78125f, 151408.21875f, 63457.26171875f, 63656.53125f, 25506.87890625f, 34634.44921875f, 42233.234375f, 44288.28515625f, 95634.6171875f, 56152.46875f,
				163879.34375f, 48217.5703125f, 191375.796875f};
		chromatograms.add(y2);
		float[] b2=new float[] {25159.392578125f, 24804.0703125f, 48990.80078125f, 17596.7421875f, 13867.8037109375f, 0.0f, 0.0f, 0.0f, 56607.82421875f, 55368.59375f, 122756.265625f, 136530.40625f,
				147656.84375f, 73210.0625f, 55667.26953125f, 20306.173828125f, 18992.7265625f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
		chromatograms.add(b2);
		float[] y3=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 60140.7890625f, 77383.78125f, 234447.40625f, 384764.125f, 551070.3125f, 365833.71875f, 340419.40625f, 185743.34375f,
				89611.265625f, 53095.55078125f, 18336.125f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
		chromatograms.add(y3);
		float[] b3=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 20709.3359375f, 0.0f, 19313.189453125f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
				0.0f};
		chromatograms.add(b3);
		float[] y4=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 127182.1640625f, 201739.125f, 245288.375f, 220938.375f, 153956.828125f, 104623.1015625f, 58068.58203125f,
				67807.0234375f, 22975.841796875f, 24826.98046875f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
		chromatograms.add(y4);
		float[] y5=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 151542.015625f, 208508.59375f, 181994.015625f, 0.0f, 153442.9375f, 149518.125f, 102380.40625f, 69725.234375f, 46574.99609375f,
				29580.294921875f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
		chromatograms.add(y5);
		float[] y6=new float[] {0.0f, 30776.84765625f, 39408.5f, 67201.9921875f, 33696.421875f, 35307.5546875f, 0.0f, 70357.0234375f, 25180.89453125f, 53548.546875f, 73943.6875f, 89488.25f,
				73816.078125f, 69158.9375f, 14915.70703125f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
		chromatograms.add(y6);
		float[] y7=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 109685.7109375f, 228564.53125f, 377972.59375f, 682637.5f, 763072.25f, 591621.0625f, 446442.8125f, 278002.375f, 150419.625f,
				69802.625f, 34787.01171875f, 19528.890625f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
		chromatograms.add(y7);
		float[] y8=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 133591.40625f, 105622.3203125f, 107714.359375f, 90343.625f, 187992.90625f, 175194.765625f,
				400795.59375f, 218788.625f, 387960.625f, 154844.1875f, 276996.125f, 0.0f, 0.0f, 0.0f};
		chromatograms.add(y8);
		String[] ionNames=new String[] {"y2", "b2", "y3", "b3", "y4", "y5", "y6", "y7", "y8",};
		float[] rts=new float[] {45.98421859741211f, 46.02557373046875f, 46.06422805786133f, 46.1038818359375f, 46.142822265625f, 46.18221664428711f, 46.2204475402832f, 46.26137161254883f,
				46.29861831665039f, 46.340152740478516f, 46.3767204284668f, 46.41716003417969f, 46.455848693847656f, 46.495277404785156f, 46.53367614746094f, 46.57311248779297f, 46.61155700683594f,
				46.6527099609375f, 46.689884185791016f, 46.73152542114258f, 46.768741607666016f, 46.81028747558594f, 46.848995208740234f, 46.88846969604492f, 46.92737579345703f};

		for (int i=0; i<chromatograms.size(); i++) {
			float[] chromatogram=SkylineSGFilter.paddedSavitzkyGolaySmooth(chromatograms.get(i));
			chromatograms.set(i, chromatogram);
		}
		
		for (int i=0; i<rts.length; i++) {
			rts[i]=rts[i]*60.0f;
		}

		TDoubleArrayList masses=new TDoubleArrayList();
		int count=0;
		for (float[] f : chromatograms) {
			masses.add(count++);
		}
		double[] fragmentMasses=masses.toArray();
		
		TransitionRefinementData data=TransitionRefiner.identifyTransitions("ASVAAQQQEEAR", fragmentMasses, chromatograms, rts, true);
		float[] correlations=data.getCorrelationArray();
		float[] integrations=data.getIntegrationArray();
		for (int i=0; i<integrations.length; i++) {
			System.out.println(ionNames[i]+"\t"+correlations[i]+"\t"+integrations[i]);
		}
		Charter.launchCharts("TITLE", TransitionRefiner.getChartPanels(data));
	}

}
