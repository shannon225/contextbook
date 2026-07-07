package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.Arrays;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.KDE;
import gnu.trove.list.array.TDoubleArrayList;

public class LocalFDR {

    public static void main(String[] args) {
		double[] pValues = { 2.21118E-76, 1.97107E-26, 2.11033E-21, 1.8278E-15, 1.7355E-10, 1.31168E-09, 3.9144E-08,
				8.5185E-08, 4.1912E-07, 4.85707E-06, 1.00519E-05, 2.86505E-05, 7.20987E-05, 0.000174663, 0.00025233,
				0.000418436, 0.000895334, 0.001086783, 0.001492685, 0.002026009, 0.002221717, 0.002857665, 0.003659834,
				0.004009534, 0.00453973, 0.005476186, 0.006072238, 0.007343078, 0.008535263, 0.009463832, 0.010663848,
				0.011383, 0.012641046, 0.013518242, 0.015080017, 0.017005261, 0.019270161, 0.019801676, 0.023097154,
				0.024773068, 0.026003488, 0.028103783, 0.029195969, 0.033288393, 0.034533787, 0.036220735, 0.03956669,
				0.042097791, 0.044056716, 0.04574696, 0.049672297, 0.051816822, 0.053759541, 0.05612858, 0.057576632,
				0.059369755, 0.063058804, 0.0648954, 0.067376553, 0.069243595, 0.070718232, 0.074554566, 0.076280209,
				0.08000529, 0.085368131, 0.087120436, 0.088372546, 0.091854129, 0.094678091, 0.096368749, 0.09826717,
				0.103328418, 0.104915159, 0.107205594, 0.110134977, 0.11068028, 0.112915341, 0.114311574, 0.11738534,
				0.118841672, 0.120803892, 0.121394733, 0.123634, 0.125563591, 0.129031896, 0.135069437, 0.138449688,
				0.140800287, 0.143240825, 0.146805724, 0.150354275, 0.154646692, 0.157429613, 0.162695018, 0.16532904,
				0.171320018, 0.174089009, 0.175399167, 0.176898847, 0.179491499, 0.181886857, 0.183837359, 0.184878923,
				0.187387941, 0.189193453, 0.191280045, 0.193319929, 0.19436764, 0.200535553, 0.20262768, 0.205906035,
				0.207244613, 0.212971643, 0.215262126, 0.217310335, 0.21909639, 0.224081816, 0.228261917, 0.229410862,
				0.234964919, 0.236800978, 0.241221103, 0.242863227, 0.242909792, 0.247892895, 0.251135298, 0.252767889,
				0.256250854, 0.262330928, 0.265919998, 0.267674573, 0.268102404, 0.271567873, 0.274960591, 0.278886495,
				0.28097437, 0.28491403, 0.28912332, 0.290740657, 0.294029141, 0.297796711, 0.301712978, 0.303533723,
				0.30673464, 0.309731257, 0.310986364, 0.313608101, 0.320946324, 0.323969788, 0.325957736, 0.328497223,
				0.334418537, 0.33541152, 0.337214449, 0.340775364, 0.344825332, 0.348799907, 0.348799907, 0.34987983,
				0.351729703, 0.357594867, 0.360284663, 0.362345746, 0.366166536, 0.374004882, 0.374718994, 0.380378145,
				0.381597752, 0.386767241, 0.38891373, 0.391010422, 0.391498452, 0.39638604, 0.399218897, 0.401612201,
				0.402817152, 0.407245286, 0.41422476, 0.417610821, 0.420144725, 0.422537417, 0.425955344, 0.428047398,
				0.429760615, 0.434056534, 0.438333226, 0.440710834, 0.442500788, 0.443451894, 0.449834096, 0.451289309,
				0.455565591, 0.461296454, 0.468123438, 0.470307299, 0.470781162, 0.475175993, 0.477622046, 0.48264483,
				0.483924805, 0.488351792, 0.488536337, 0.496552294, 0.49948866, 0.500561527, 0.504294122, 0.505081347,
				0.511301717, 0.511335116, 0.511729881, 0.517640593, 0.522196878, 0.523003113, 0.525655519, 0.526319318,
				0.532629671, 0.535706909, 0.53926376, 0.543101757, 0.548261543, 0.551170067, 0.551810485, 0.559977329,
				0.561904639, 0.565204238, 0.56808267, 0.572337793, 0.575757219, 0.579090793, 0.582292121, 0.584053487,
				0.588534311, 0.589956591, 0.594139427, 0.596556454, 0.602259136, 0.603021407, 0.607530646, 0.609910957,
				0.614549601, 0.617138035, 0.617656404, 0.620604291, 0.623860076, 0.625297315, 0.629176672, 0.630344981,
				0.635194166, 0.635760896, 0.636201843, 0.639594439, 0.643899842, 0.643899842, 0.648138258, 0.650777565,
				0.651516759, 0.652843341, 0.658665105, 0.665390982, 0.671678752, 0.677730219, 0.681607557, 0.683409071,
				0.684722184, 0.688797902, 0.691095892, 0.69354455, 0.693920567, 0.697865686, 0.699577613, 0.706733145,
				0.707892405, 0.709172119, 0.711551185, 0.712299187, 0.716139103, 0.719638003, 0.720514637, 0.72469367,
				0.727411863, 0.728285726, 0.729911631, 0.732517901, 0.735830504, 0.737454203, 0.737966602, 0.739728896,
				0.742779566, 0.747932999, 0.75192639, 0.755175999, 0.757079693, 0.760737607, 0.764864879, 0.765936446,
				0.767650377, 0.770810305, 0.77488387, 0.775391682, 0.779708704, 0.783785846, 0.783785846, 0.785838679,
				0.787644529, 0.79087825, 0.791303586, 0.797074213, 0.798087485, 0.799820902, 0.805181393, 0.807719253,
				0.809879538, 0.815049633, 0.818043484, 0.819819273, 0.822069618, 0.828345744, 0.834842597, 0.837473926,
				0.839130625, 0.842297305, 0.84306888, 0.850525888, 0.856728365, 0.860019959, 0.863316755, 0.865273548,
				0.867907439, 0.87138183, 0.871512837, 0.874595598, 0.874710089, 0.877579771, 0.877806688, 0.880360441,
				0.882959727, 0.885396525, 0.887687085, 0.889845467, 0.891883913, 0.893432674, 0.895642534, 0.898227323,
				0.901369223, 0.902115082, 0.905248204, 0.90779727, 0.909866041, 0.912108172, 0.918725602, 0.923345405,
				0.924718434, 0.925377796, 0.929553316, 0.931645149, 0.936532286, 0.940539827, 0.946668305, 0.956276285,
				0.999934206, 0.999956091, 0.99996371, 0.999967899, 0.999971237, 0.999972383, 0.999976147, 0.999977541,
				0.999978378, 0.999979471, 0.999980156, 0.999981602, 0.99998211, 0.999982633, 0.9999829, 0.999983448,
				0.999984305, 0.999984601, 0.999984601, 0.999985211, 0.999985526, 0.999985848, 0.999986177, 0.999986514,
				0.999986514, 0.99998686, 0.999987215, 0.999987581, 0.999988346, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        double[] lfdrValues = estimateLocalFDR(pValues);
        
        double[] fdrValues=BenjaminiHochberg.calculateAdjustedPValues(pValues);

        for (int i = 0; i < pValues.length; i++) {
            System.out.println(pValues[i]+"\t"+fdrValues[i]+"\t"+lfdrValues[i]);
		}
    }

    /**
     * based on Logit from https://github.com/StoreyLab/qvalue/blob/master/R/lfdr.R
     * @param pvalues
     * @return
     */
    public static double[] estimateLocalFDR(double[] pvalues) {
    	
        double[] lfdr = new double[pvalues.length];
        double smallestMeaningfulP=1e-8; // estimate from Storey et al

        double pi0 = estimatePi0(pvalues);

        int n = pvalues.length;

        // logit transform
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = Math.log((pvalues[i] + smallestMeaningfulP) / (1 - pvalues[i] + smallestMeaningfulP));
        }
        KDE density;
        try {
        	density=new KDE(x, 1.0, 100);
        } catch (Exception e) {
        	Logger.errorLine("ERRORED WITH "+pvalues.length+" VALUES: "+General.toString(pvalues));
        	density=null; // FIXME
        	throw e;
        }
        
        for (int i = 0; i < n; i++) {
            double dx = Math.exp(x[i]) / Math.pow(1 + Math.exp(x[i]), 2);
            lfdr[i] = (pi0 * dx) / density.getProbability(x[i]);
        }

        // truncate
        for (int i = 0; i < lfdr.length; i++) {
            lfdr[i] = Math.min(lfdr[i], 1);
        }

        // make monotonic increasing
        Integer[] order = new Integer[n];
        for (int i = 0; i < n; i++) {
            order[i] = i;
        }
        Arrays.sort(order, (a, b) -> Double.compare(pvalues[a], pvalues[b]));

        double[] lfdrOrdered = new double[n];
        for (int i = 0; i < n; i++) {
            lfdrOrdered[i] = lfdr[order[i]];
        }

        for (int i = 1; i < n; i++) {
            lfdrOrdered[i] = Math.max(lfdrOrdered[i], lfdrOrdered[i - 1]);
        }

        for (int i = 0; i < n; i++) {
            lfdr[order[i]] = lfdrOrdered[i];
        }

        return lfdr;
    }

    /**
     * Storey-style estimate of the proportion of true null hypotheses.
     * Returns the median of per-lambda estimates, then caps at 1.0 (the
     * theoretical upper bound). Without the cap the estimator can exceed 1
     * when sample sizes are small or the p-value distribution is noisy in
     * the upper range -- values > 1 are mathematically meaningless and
     * propagate into q-values / local FDRs > 1 downstream.
     */
    public static double estimatePi0(double[] pValues) {
    	double[] lambdas=new double[9];
    	double[] pi0s=new double[lambdas.length];
    	for (int i = 0; i < lambdas.length; i++) {
			lambdas[i]=(i+1)/(double)(lambdas.length+1);
			pi0s[i]=estimatePi0Local(pValues, lambdas[i]);
		}

    	return Math.min(1.0, QuickMedianDouble.median(pi0s));
    }

    /**
     * Per-lambda pi_0 estimate. Not capped -- callers using this directly
     * should cap themselves. Public callers should prefer {@link #estimatePi0}.
     */
    public static double estimatePi0Local(double[] pValues, double lambda) {
        int m = pValues.length;
        int count = 0;
        for (double p : pValues) {
            if (p > lambda) {
                count++;
            }
        }
        return (double) count / (m * (1 - lambda));
    }

}
