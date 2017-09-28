package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;

public class TransitionRefinerTest {

	public static void main(String[] args) throws Exception {
		//plotEIGNIISDAMK();
		
		String sequence="AAHSAELEAVLLALAR";
		byte precursorCharge=(byte)3;
		float targetRT=87.86211f*60f;

		File diaFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/zero_hela/121115_bcs_hela_24mz_400_1000_0D_1.dia");
		HashMap<String, String> defaults=SearchParameterParser.getDefaultParameters();
		SearchParameters params=SearchParameterParser.parseParameters(defaults);
		float rtRange=params.getExpectedPeakWidth();
		StripeFileInterface dia=StripeFileGenerator.getFile(diaFile, params, true);

		FragmentationModel model=new FragmentationModel(sequence, params.getAAConstants());
		AnnotatedLibraryEntry entry=model.getUnitSpectrum(diaFile.getName(), new HashSet<>(), precursorCharge, targetRT, params);
		ArrayList<Stripe> stripes=dia.getStripes(entry.getPrecursorMZ(), targetRT-rtRange, targetRT+rtRange, false);
		Collections.sort(stripes);

		PSMData psmdata=new PSMData(entry.getAccessions(), entry.getSpectrumIndex(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), entry.getPeptideModSeq(), targetRT, entry.getScore(), 1.0f-entry.getScore(), 2*rtRange, false);
		PeptideQuantExtractorTask quantTask=new PeptideQuantExtractorTask(dia.getOriginalFileName(), psmdata, Optional.empty(), stripes, params, false);
		TransitionRefinementData data=quantTask.extractSpectrum(entry, rtRange, false, false);
		Charter.launchCharts("TITLE", TransitionRefiner.getChartPanels(data));
	}
	
	//EIGNIISDAMK+2 (rt=4229.2856)
	public static void plotEIGNIISDAMK() {
		ArrayList<float[]> chromatograms=new ArrayList<float[]>();
		float[] b2=new float[] {83875.0390625f, 102077.375f, 56852.09765625f, 71273.8515625f, 56333.8359375f, 72436.3984375f, 82192.1640625f, 209511.0f, 186783.5f, 294653.34375f, 363570.65625f,
				272646.5625f, 288544.3125f, 292669.5f, 183966.25f, 166379.015625f, 162588.015625f, 137269.28125f, 141212.984375f, 259367.15625f, 348729.46875f, 114987.4140625f, 81208.6953125f,
				67690.3125f, 84733.7578125f, 112416.140625f,};
		chromatograms.add(b2);
		float[] y2=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 102044.9921875f, 183339.03125f, 186134.890625f, 269246.34375f, 211908.0625f, 125531.7890625f, 0.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 0.0f, 3050.194091796875f, 25618.966796875f, 27664.83203125f, 34760.046875f, 21478.712890625f,};
		chromatograms.add(y2);
		float[] b3=new float[] {28463.4921875f, 30860.783203125f, 25471.0859375f, 35515.40234375f, 36057.078125f, 28891.68359375f, 30550.974609375f, 41254.03125f, 32727.703125f, 41228.16796875f, 0.0f,
				30068.244140625f, 49734.046875f, 75036.4296875f, 132127.171875f, 247096.96875f, 345049.8125f, 492589.625f, 350056.84375f, 423721.15625f, 424324.03125f, 155052.625f, 59360.97265625f,
				53621.98046875f, 44607.37890625f, 65869.4296875f,};
		chromatograms.add(b3);
		float[] y3=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 88795.4765625f, 335280.75f, 297209.59375f, 386607.625f, 335113.21875f, 204680.6875f, 105236.3359375f,
				41092.90625f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,};
		chromatograms.add(y3);
		float[] b4=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 76503.7734375f, 124448.03125f, 84324.5078125f, 47245.15625f, 59678.08984375f, 88522.765625f,
				58771.078125f, 144687.5f, 32791.0078125f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,};
		chromatograms.add(b4);
		float[] y4=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 99901.8203125f, 158237.609375f, 145698.984375f, 86498.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
				12157.671875f, 0.0f, 0.0f, 0.0f, 0.0f,};
		chromatograms.add(y4);
		float[] b5=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 35020.3515625f, 42871.60546875f, 30500.185546875f,
				0.0f, 0.0f, 0.0f, 0.0f,};
		chromatograms.add(b5);
		float[] y5=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 262998.5f, 687530.25f, 700176.0625f, 1041231.75f, 747892.75f, 485779.3125f, 212630.515625f, 142454.09375f,
				45828.2890625f, 39074.15625f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,};
		chromatograms.add(y5);
		float[] y6=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 13842.5771484375f, 41806.62890625f, 0.0f, 300604.78125f, 664254.625f, 673934.0625f, 898653.125f, 793304.0f, 514485.0625f,
				358957.46875f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,};
		chromatograms.add(y6);
		float[] y7=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 19504.962890625f, 83693.8203125f, 147224.3125f, 157680.515625f, 281689.28125f, 247816.375f, 144363.15625f, 0.0f, 0.0f,
				0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,};
		chromatograms.add(y7);
		float[] y8=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 27307.0f, 127635.953125f, 45227.5625f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
				0.0f, 0.0f,};
		chromatograms.add(y8);
		float[] y9=new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 168510.5625f, 393503.1875f, 368990.21875f, 611638.8125f, 419877.75f, 324253.5625f, 173974.625f, 90743.15625f,
				0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 8127.353515625f, 5719.25244140625f, 0.0f, 0.0f,};
		chromatograms.add(y9);
		float[] rts=new float[] {70.01758575439453f, 70.0562744140625f, 70.09734344482422f, 70.1347885131836f, 70.17521667480469f, 70.21282196044922f, 70.25397491455078f, 70.293701171875f,
				70.3323974609375f, 70.37088775634766f, 70.40955352783203f, 70.44808959960938f, 70.48809051513672f, 70.52470397949219f, 70.5665054321289f, 70.60203552246094f, 70.64164733886719f,
				70.67855072021484f, 70.71800994873047f, 70.7557373046875f, 70.79398345947266f, 70.83085632324219f, 70.87123107910156f, 70.9073257446289f, 70.94801330566406f, 70.98357391357422f,};

		/*for (int i=0; i<chromatograms.size(); i++) {
			float[] chromatogram=SkylineSGFilter.paddedSavitzkyGolaySmooth(chromatograms.get(i));
			chromatograms.set(i, chromatogram);
		}*/

		for (int i=0; i<rts.length; i++) {
			rts[i]=rts[i]*60.0f;
		}

		String[] ionNames=new String[] {"b2", "y2", "b3", "y3", "b4", "y4", "b5", "y5", "y6", "y7", "y8", "y9",};
		ArrayList<FragmentIon> ions=new ArrayList<FragmentIon>();
		ions.add(new FragmentIon(2, (byte)2, IonType.b));
		ions.add(new FragmentIon(2, (byte)2, IonType.y));
		ions.add(new FragmentIon(3, (byte)3, IonType.b));
		ions.add(new FragmentIon(3, (byte)3, IonType.y));
		ions.add(new FragmentIon(4, (byte)4, IonType.b));
		ions.add(new FragmentIon(4, (byte)4, IonType.y));
		ions.add(new FragmentIon(5, (byte)5, IonType.b));
		ions.add(new FragmentIon(5, (byte)5, IonType.y));
		ions.add(new FragmentIon(6, (byte)6, IonType.y));
		ions.add(new FragmentIon(7, (byte)7, IonType.y));
		ions.add(new FragmentIon(8, (byte)8, IonType.y));
		ions.add(new FragmentIon(9, (byte)9, IonType.y));
		FragmentIon[] fragmentMasses=ions.toArray(new FragmentIon[ions.size()]);

		TransitionRefinementData data=TransitionRefiner.identifyTransitions("EIGNIISDAMK", (byte)2, 70.40955352783203f, fragmentMasses, chromatograms, rts, Optional.ofNullable((float[])null), true);
		float[] correlations=data.getCorrelationArray();
		float[] integrations=data.getIntegrationArray();
		for (int i=0; i<integrations.length; i++) {
			System.out.println(ionNames[i]+"\t"+correlations[i]+"\t"+integrations[i]);
		}
		Charter.launchCharts("TITLE", TransitionRefiner.getChartPanels(data));
	}

	public static void plotASVAAQQQEEAR() {
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

		String[] ionNames=new String[] {"y2", "b2", "y3", "b3", "y4", "y5", "y6", "y7", "y8",};
		ArrayList<FragmentIon> ions=new ArrayList<FragmentIon>();
		ions.add(new FragmentIon(2, (byte)2, IonType.y));
		ions.add(new FragmentIon(2, (byte)2, IonType.b));
		ions.add(new FragmentIon(3, (byte)3, IonType.y));
		ions.add(new FragmentIon(3, (byte)3, IonType.b));
		ions.add(new FragmentIon(4, (byte)4, IonType.y));
		ions.add(new FragmentIon(5, (byte)5, IonType.y));
		ions.add(new FragmentIon(6, (byte)6, IonType.y));
		ions.add(new FragmentIon(7, (byte)7, IonType.y));
		ions.add(new FragmentIon(8, (byte)8, IonType.y));
		FragmentIon[] fragmentMasses=ions.toArray(new FragmentIon[ions.size()]);

		TransitionRefinementData data=TransitionRefiner.identifyTransitions("ASVAAQQQEEAR", (byte)2, 46.41716003417969f, fragmentMasses, chromatograms, rts, Optional.ofNullable((float[])null), true);
		float[] correlations=data.getCorrelationArray();
		float[] integrations=data.getIntegrationArray();
		for (int i=0; i<integrations.length; i++) {
			System.out.println(ionNames[i]+"\t"+correlations[i]+"\t"+integrations[i]);
		}
		Charter.launchCharts("TITLE", TransitionRefiner.getChartPanels(data));
	}

}
