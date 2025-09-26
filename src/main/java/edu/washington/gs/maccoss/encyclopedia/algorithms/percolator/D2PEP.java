package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.MonotonicCubicSplineInterpolatedFunction;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;

class ScoredPSM implements Comparable<ScoredPSM> {
    private final double primaryScore;
    private final double estimatedPEP;
    private final boolean decoy;
    public ScoredPSM(double primaryScore, boolean decoy, double estimatedPEP) {
        this.primaryScore = primaryScore;
        this.decoy = decoy;
        this.estimatedPEP=estimatedPEP;
    }
    public double getPrimaryScore() { return primaryScore; }
    public boolean isDecoy() { return decoy; }
    public double getEstimatedPEP() { return estimatedPEP; }
    
    @Override
    public int compareTo(ScoredPSM o) {
		if (o==null) return 1;
		int c=Double.compare(primaryScore, o.primaryScore);
		if (c!=0) return c;
		if (decoy==o.decoy) return 0;
		if (decoy) return 1;
		return -1;
    }
}

public class D2PEP {
    private static final boolean ASSUME_TDC = false;// with equal targets and decoys
    
    // Configurable parameters for adaptive binning
    public static final int DEFAULT_NUM_BINS = 1000;
    public static final int DEFAULT_MIN_DECOYS = 3;
    public static final int DEFAULT_MAX_BIN_WIDTH = 20;

	static class ScoredLabel {
        double score;
        int label; // 0=target, 1=decoy

        ScoredLabel(double score, int label) {
            this.score = score;
            this.label = label;
        }
    }

    /**
     * Main function: returns PEPs for all observations (in input order).
     * Based on Yuqi Zheng's code:
     * https://github.com/statisticalbiotechnology/smooth_q_to_pep/blob/main/pyIsoPEP/IsotonicPEP.py
     */
    public static Function d2pep(List<ScoredPSM> psms) {
        int n = psms.size();

        // Step 1: Collect ObsInfo, sort by score descending
        List<ScoredLabel> obsList = new ArrayList<>();
        double maxScore=Float.NEGATIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            ScoredPSM psm = psms.get(i);
            obsList.add(new ScoredLabel(psm.getPrimaryScore(), psm.isDecoy() ? 1 : 0));
            maxScore=Math.max(maxScore, psm.getPrimaryScore());
        }
        obsList.sort(Comparator.comparingDouble((ScoredLabel o) -> o.score).reversed());

        // Step 2: Add pseudo-observation at the top
        ScoredLabel pseudo = new ScoredLabel(maxScore+Math.abs(maxScore*0.01), 0); // Label=0.5 is handled below 1% higher than max score
        List<ScoredLabel> withPseudo = new ArrayList<>();
        withPseudo.add(pseudo);
        withPseudo.addAll(obsList);

        // Step 3: Prepare arrays for binning (score, label)
        double[] scores = new double[withPseudo.size()];
        double[] labels = new double[withPseudo.size()];
        for (int i = 0; i < withPseudo.size(); i++) {
            scores[i] = withPseudo.get(i).score;
            labels[i] = (i == 0) ? 0.5 : withPseudo.get(i).label; // pseudo=0.5, else 0/1
        }

        // Step 4: Adaptive binning (returns bin centers, avg label, weight)
        BinnedData binned = adaptiveBinning(scores, labels, DEFAULT_NUM_BINS, DEFAULT_MAX_BIN_WIDTH, DEFAULT_MIN_DECOYS);

        // --- Kernel smoothing ---
        double medianDx = 0.0;
        for (int i = 1; i < binned.x.length; i++) medianDx += Math.abs(binned.x[i] - binned.x[i-1]);
        medianDx /= (binned.x.length - 1);
        double bandwidth = medianDx * 2.0; // or tune empirically

        double[] smoothY = kernelSmooth(binned.x, binned.y, bandwidth);

        // Step 5: Storey Pool Adjacent Violators Algorithm (PAVA)
        double[] pepFit = isotonicRegression(binned.x, smoothY);
        
        double pi0;
        if (ASSUME_TDC) {
            // currently assumes #T=#D for T/D competition!
        	pi0=1.0;
        } else {
            double worstP=General.max(pepFit);
            pi0=(1-worstP)/worstP;
        }
		//System.out.println("pi0: "+pi0);
		
        ArrayList<XYPoint> knots=new ArrayList<XYPoint>();
        for (int i = 0; i < binned.x.length; i++) {
        	double p=Math.min(1.0, pi0*pepFit[i] / (1.0 - pepFit[i]));
        	knots.add(new XYPoint(binned.x[i], p));
		}
        Collections.sort(knots);

        MonotonicCubicSplineInterpolatedFunction function=new MonotonicCubicSplineInterpolatedFunction(knots, false);
        
//        for (int i = 0; i < binned.x.length; i++) {
//			System.out.println(binned.x[i]+"\t"+binned.y[i]+"\t"+pepFit[i]+"\t"+function.getYValue((float)binned.x[i]));
//		}
        
        return function;
    }

    // Holds binned representation: x (score position), y (mean label), w (count)
    static class BinnedData {
        double[] x; double[] y; int[] w;
        BinnedData(double[] x, double[] y, int[] w) { this.x = x; this.y = y; this.w = w; }
    }

    static BinnedData adaptiveBinning(
            double[] scores, double[] labels, int maxBins, int binWidth, int unusedMinDecoys) {

        int n = scores.length;
        int numBins = Math.max(1, (int)Math.round((double)n / binWidth));
        numBins = Math.min(numBins, maxBins);

        int minBinSize = n / numBins; // floor division
        int remainder = n % numBins;  // Number of bins that get +1

        
        TDoubleArrayList bx = new TDoubleArrayList();
        TDoubleArrayList by = new TDoubleArrayList();
        TIntArrayList bw = new TIntArrayList();

        int start = 0;
        for (int b = 0; b < numBins; b++) {
            int thisBinSize = minBinSize + (b < remainder ? 1 : 0);
            int end = start + thisBinSize;

            double meanPos = mean(scores, start, end);
            double meanLab = mean(labels, start, end);
            bx.add(meanPos);
            by.add(meanLab);
            bw.add(end - start);

            start = end;
        }
        double[] bxArr = bx.toArray();
        double[] byArr = by.toArray();
        int[] bwArr = bw.toArray();
        return new BinnedData(bxArr, byArr, bwArr);
    }

    // Helper: mean of array slice [start, end)
    static double mean(double[] arr, int start, int end) {
        double s = 0;
        int cnt = end - start;
        for (int i = start; i < end; i++) s += arr[i];
        return s / cnt;
    }
    
    static double[] kernelSmooth(double[] x, double[] y, double bandwidth) {
        int n = x.length;
        double[] ySmooth = new double[n];
        double denom = 2 * bandwidth * bandwidth;
        for (int i = 0; i < n; i++) {
            double num = 0.0;
            double weightSum = 0.0;
            for (int j = 0; j < n; j++) {
                double dist = x[i] - x[j];
                double w = Math.exp(-dist * dist / denom);
                num += y[j] * w;
                weightSum += w;
            }
            ySmooth[i] = num / weightSum;
        }
        return ySmooth;
    }
    
    public static double[] isotonicRegression(double[] x, double[] y) {
        int n = y.length;
        double[] fitted = y.clone();
        int[] block = new int[n];
        double[] sum = new double[n];
        int[] count = new int[n];
        for (int i = 0; i < n; i++) {
            block[i] = i;
            sum[i] = y[i];
            count[i] = 1;
        }
        for (int i = 0; i < n - 1; ) {
            if (fitted[i] > fitted[i+1]) {
                int j = i;
                while (j >= 0 && fitted[j] > fitted[j+1]) {
                    // Merge blocks
                    int left = block[j];
                    int right = block[j+1];
                    double mergedSum = sum[left] + sum[right];
                    int mergedCount = count[left] + count[right];
                    double mergedVal = mergedSum / mergedCount;
                    for (int k = left; k <= right; k++) {
                        fitted[k] = mergedVal;
                        block[k] = left;
                        sum[left] = mergedSum;
                        count[left] = mergedCount;
                    }
                    j--;
                }
            } else {
                i++;
            }
        }

        for (int i = 0; i < fitted.length-1; i++) {
            if (fitted[i] > fitted[i+1]) {
            	fitted[i+1]=fitted[i];
            }
        }
        return fitted;
    }
    
    /**
     * I-spline style monotonic interpolation using a cumulative smooth-step basis.
     * Matches the Python behavior:
     *
     *   f(x) = yB[0] + sum_{j=1}^{n-1} (yB[j] - yB[j-1]) * I_j(x)
     *
     * where
     *   I_j(x) = 0                       if x <= xB[j-1]
     *            3u^2 - 2u^3             if xB[j-1] < x < xB[j]
     *            1                       if x >= xB[j]
     *
     * with u = (x - xB[j-1]) / (xB[j] - xB[j-1]).
     *
     * Requirements:
     *  - xB must be strictly increasing (ties are skipped safely).
     *  - yB should be monotone (non-decreasing for increasing f, or non-increasing for decreasing f).
     *
     * @param xB    knot (block-center) x positions, length n >= 0
     * @param yB    knot values (block means), length n
     * @param xEval evaluation points, length m
     * @return interpolated values f(xEval), length m
     */
    public static double[] iSplineMonotonicInterpolate(double[] xB, double[] yB, double[] xEval) {
        final int n = xB == null ? 0 : xB.length;

        if (n == 0) {
            return new double[0]; // mirrors the Python function
        }
        if (n == 1) {
            double[] out = new double[xEval.length];
            Arrays.fill(out, yB[0]);
            return out;
        }

        // Defensive copy not strictly needed; do it if you might mutate upstream arrays.
        xB = Arrays.copyOf(xB, n);
        yB = Arrays.copyOf(yB, n);

        // Precompute segment deltas (allowing positive OR negative for either monotone direction)
        final double[] d = new double[n - 1];
        for (int j = 1; j < n; j++) {
            d[j - 1] = yB[j] - yB[j - 1];
        }

        final int m = xEval.length;
        final double[] f = new double[m];
        Arrays.fill(f, yB[0]); // base level

        // Accumulate smooth cumulative contributions for each interval
        for (int j = 1; j < n; j++) {
            final double xl = xB[j - 1];
            final double xr = xB[j];
            final double dx = xr - xl;

            if (!(dx > 0)) {
                // skip zero-width or invalid intervals
                continue;
            }

            final double dj = d[j - 1];

            // For each xEval, add dj * I_j(x)
            // I_j(x) = smoothstep(u) with u in [0,1], clipped outside.
            for (int i = 0; i < m; i++) {
                double u = (xEval[i] - xl) / dx;
                if (u <= 0.0) {
                    // I_j = 0
                    continue;
                } else if (u >= 1.0) {
                    // I_j = 1
                    f[i] += dj;
                } else {
                    // cubic smoothstep: 3u^2 - 2u^3
                    double u2 = u * u;
                    double s = (3.0 * u2) - (2.0 * u2 * u);
                    f[i] += dj * s;
                }
            }
        }

        return f;
    }

    // Example usage
    public static void main(String[] args) {
        List<ScoredPSM> psms = Arrays.asList(
        		new ScoredPSM(3.42083,false,0.00000117599),
        		new ScoredPSM(3.39605,false,0.00000128611),
        		new ScoredPSM(3.36827,false,0.00000142182),
        		new ScoredPSM(3.24497,false,0.00000221955),
        		new ScoredPSM(3.193,false,0.00000267788),
        		new ScoredPSM(3.17951,false,0.00000281163),
        		new ScoredPSM(3.05944,false,0.00000433818),
        		new ScoredPSM(2.9489,false,0.0000064671),
        		new ScoredPSM(2.88734,false,0.00000807732),
        		new ScoredPSM(2.87771,false,0.00000836324),
        		new ScoredPSM(2.83641,false,0.00000970864),
        		new ScoredPSM(2.819,false,0.0000103391),
        		new ScoredPSM(2.81245,false,0.0000105865),
        		new ScoredPSM(2.77004,false,0.0000123388),
        		new ScoredPSM(2.75336,false,0.0000131051),
        		new ScoredPSM(2.74951,false,0.0000132886),
        		new ScoredPSM(2.74019,false,0.0000137438),
        		new ScoredPSM(2.70989,false,0.000015333),
        		new ScoredPSM(2.70956,false,0.0000153518),
        		new ScoredPSM(2.69797,false,0.0000160079),
        		new ScoredPSM(2.6872,false,0.0000166429),
        		new ScoredPSM(2.67692,false,0.0000172721),
        		new ScoredPSM(2.64509,false,0.0000193769),
        		new ScoredPSM(2.61966,false,0.0000212407),
        		new ScoredPSM(2.6169,false,0.0000214534),
        		new ScoredPSM(2.60716,false,0.0000222223),
        		new ScoredPSM(2.6064,false,0.0000222828),
        		new ScoredPSM(2.59526,false,0.000023198),
        		new ScoredPSM(2.59295,false,0.0000233919),
        		new ScoredPSM(2.59032,false,0.0000236155),
        		new ScoredPSM(2.55681,false,0.0000266541),
        		new ScoredPSM(2.55644,false,0.0000266901),
        		new ScoredPSM(2.55111,false,0.0000272086),
        		new ScoredPSM(2.54963,false,0.0000273543),
        		new ScoredPSM(2.53915,false,0.0000284101),
        		new ScoredPSM(2.53777,false,0.0000285517),
        		new ScoredPSM(2.53636,false,0.0000286974),
        		new ScoredPSM(2.53622,false,0.0000287119),
        		new ScoredPSM(2.53552,false,0.0000287846),
        		new ScoredPSM(2.53409,false,0.0000289338),
        		new ScoredPSM(2.53276,false,0.0000290732),
        		new ScoredPSM(2.53239,false,0.0000291116),
        		new ScoredPSM(2.49528,false,0.000033288),
        		new ScoredPSM(2.48407,false,0.0000346636),
        		new ScoredPSM(2.48289,false,0.000034811),
        		new ScoredPSM(2.4766,false,0.0000356117),
        		new ScoredPSM(2.47319,false,0.0000360527),
        		new ScoredPSM(2.47308,false,0.0000360668),
        		new ScoredPSM(2.47275,false,0.0000361099),
        		new ScoredPSM(2.4653,false,0.0000370947),
        		new ScoredPSM(2.45971,false,0.0000378519),
        		new ScoredPSM(2.45248,false,0.0000388532),
        		new ScoredPSM(2.45234,false,0.0000388722),
        		new ScoredPSM(2.42933,false,0.000042241),
        		new ScoredPSM(2.4291,false,0.0000422768),
        		new ScoredPSM(2.42827,false,0.0000424033),
        		new ScoredPSM(2.42416,false,0.0000430381),
        		new ScoredPSM(2.42209,false,0.0000433607),
        		new ScoredPSM(2.41465,false,0.0000445417),
        		new ScoredPSM(2.40707,false,0.0000457784),
        		new ScoredPSM(2.40311,false,0.0000464378),
        		new ScoredPSM(2.39957,false,0.0000470349),
        		new ScoredPSM(2.38033,false,0.0000504197),
        		new ScoredPSM(2.37823,false,0.0000508038),
        		new ScoredPSM(2.37769,false,0.0000509027),
        		new ScoredPSM(2.37677,false,0.0000510731),
        		new ScoredPSM(2.37591,false,0.0000512305),
        		new ScoredPSM(2.3671,false,0.0000528878),
        		new ScoredPSM(2.35708,false,0.0000548367),
        		new ScoredPSM(2.35443,false,0.000055365),
        		new ScoredPSM(2.35246,false,0.0000557601),
        		new ScoredPSM(2.34497,false,0.0000572888),
        		new ScoredPSM(2.33888,false,0.0000585639),
        		new ScoredPSM(2.33221,false,0.0000599899),
        		new ScoredPSM(2.32183,false,0.0000622821),
        		new ScoredPSM(2.31906,false,0.0000629101),
        		new ScoredPSM(2.30456,false,0.0000662917),
        		new ScoredPSM(2.30233,false,0.0000668286),
        		new ScoredPSM(2.30117,false,0.0000671092),
        		new ScoredPSM(2.28969,false,0.0000699494),
        		new ScoredPSM(2.28277,false,0.0000717199),
        		new ScoredPSM(2.27666,false,0.0000733195),
        		new ScoredPSM(2.26686,false,0.0000759616),
        		new ScoredPSM(2.26227,false,0.0000772316),
        		new ScoredPSM(2.26043,false,0.0000777474),
        		new ScoredPSM(2.23304,false,0.000085831),
        		new ScoredPSM(2.21971,false,0.0000900637),
        		new ScoredPSM(2.2178,false,0.0000906881),
        		new ScoredPSM(2.21059,false,0.000093083),
        		new ScoredPSM(2.20992,false,0.0000933079),
        		new ScoredPSM(2.20375,false,0.00009541),
        		new ScoredPSM(2.20302,false,0.0000956607),
        		new ScoredPSM(2.19513,false,0.0000984267),
        		new ScoredPSM(2.18572,false,0.000101829),
        		new ScoredPSM(2.1763,false,0.000105353),
        		new ScoredPSM(2.17161,false,0.000107155),
        		new ScoredPSM(2.16039,false,0.000111587),
        		new ScoredPSM(2.14833,false,0.000116554),
        		new ScoredPSM(2.12448,false,0.000127039),
        		new ScoredPSM(2.11814,false,0.000129983),
        		new ScoredPSM(2.11521,false,0.000131365),
        		new ScoredPSM(2.11219,false,0.000132803),
        		new ScoredPSM(2.11112,false,0.000133318),
        		new ScoredPSM(2.09866,false,0.000139456),
        		new ScoredPSM(2.0947,false,0.000141463),
        		new ScoredPSM(2.09106,false,0.000143338),
        		new ScoredPSM(2.08517,false,0.000146417),
        		new ScoredPSM(2.07761,false,0.00015047),
        		new ScoredPSM(2.07747,false,0.000150549),
        		new ScoredPSM(2.06047,false,0.000160081),
        		new ScoredPSM(2.05196,false,0.000165074),
        		new ScoredPSM(2.04806,false,0.00016742),
        		new ScoredPSM(2.04625,false,0.000168516),
        		new ScoredPSM(2.04596,false,0.000168695),
        		new ScoredPSM(2.03359,false,0.000176401),
        		new ScoredPSM(2.00799,false,0.000193488),
        		new ScoredPSM(2.00025,false,0.000198975),
        		new ScoredPSM(1.98575,false,0.00020967),
        		new ScoredPSM(1.98269,false,0.000211997),
        		new ScoredPSM(1.97724,false,0.000216213),
        		new ScoredPSM(1.97503,false,0.00021795),
        		new ScoredPSM(1.97272,false,0.000219771),
        		new ScoredPSM(1.96632,false,0.000224907),
        		new ScoredPSM(1.96559,false,0.000225507),
        		new ScoredPSM(1.95989,false,0.000230197),
        		new ScoredPSM(1.95815,false,0.000231643),
        		new ScoredPSM(1.95068,false,0.000237981),
        		new ScoredPSM(1.93789,false,0.000249228),
        		new ScoredPSM(1.93352,false,0.000253198),
        		new ScoredPSM(1.93342,false,0.000253283),
        		new ScoredPSM(1.93155,false,0.000254999),
        		new ScoredPSM(1.92403,false,0.000262025),
        		new ScoredPSM(1.92328,false,0.00026273),
        		new ScoredPSM(1.9225,false,0.000263479),
        		new ScoredPSM(1.92009,false,0.000265777),
        		new ScoredPSM(1.91505,false,0.000270656),
        		new ScoredPSM(1.90959,false,0.000276045),
        		new ScoredPSM(1.90868,false,0.00027696),
        		new ScoredPSM(1.90765,false,0.000277993),
        		new ScoredPSM(1.90628,false,0.000279365),
        		new ScoredPSM(1.90144,false,0.000284294),
        		new ScoredPSM(1.90132,false,0.000284418),
        		new ScoredPSM(1.89769,false,0.000288171),
        		new ScoredPSM(1.89665,false,0.000289253),
        		new ScoredPSM(1.89282,false,0.000293283),
        		new ScoredPSM(1.8928,false,0.000293302),
        		new ScoredPSM(1.8924,false,0.000293728),
        		new ScoredPSM(1.89014,false,0.000296133),
        		new ScoredPSM(1.88773,false,0.000298729),
        		new ScoredPSM(1.8828,false,0.00030409),
        		new ScoredPSM(1.87988,false,0.000307315),
        		new ScoredPSM(1.87534,false,0.000312399),
        		new ScoredPSM(1.87454,false,0.000313301),
        		new ScoredPSM(1.86742,false,0.000321458),
        		new ScoredPSM(1.86486,false,0.000324445),
        		new ScoredPSM(1.85814,false,0.00033241),
        		new ScoredPSM(1.85355,false,0.00033797),
        		new ScoredPSM(1.84074,false,0.000353973),
        		new ScoredPSM(1.83491,false,0.000361503),
        		new ScoredPSM(1.83248,false,0.000364692),
        		new ScoredPSM(1.83021,false,0.000367687),
        		new ScoredPSM(1.82208,false,0.000378637),
        		new ScoredPSM(1.82176,false,0.000379078),
        		new ScoredPSM(1.81346,false,0.000390613),
        		new ScoredPSM(1.81317,false,0.000391026),
        		new ScoredPSM(1.80446,false,0.000403518),
        		new ScoredPSM(1.78482,false,0.00043317),
        		new ScoredPSM(1.75125,false,0.000488992),
        		new ScoredPSM(1.75044,false,0.000490422),
        		new ScoredPSM(1.74852,false,0.000493835),
        		new ScoredPSM(1.74134,false,0.000506803),
        		new ScoredPSM(1.73864,false,0.000511773),
        		new ScoredPSM(1.72798,false,0.000531843),
        		new ScoredPSM(1.72163,false,0.000544179),
        		new ScoredPSM(1.72124,false,0.00054495),
        		new ScoredPSM(1.7188,false,0.000549778),
        		new ScoredPSM(1.71306,false,0.000561282),
        		new ScoredPSM(1.70156,false,0.000585074),
        		new ScoredPSM(1.69009,false,0.000609802),
        		new ScoredPSM(1.68974,false,0.000610587),
        		new ScoredPSM(1.68128,false,0.000629511),
        		new ScoredPSM(1.66429,false,0.000669338),
        		new ScoredPSM(1.65915,false,0.000681877),
        		new ScoredPSM(1.65566,false,0.000690508),
        		new ScoredPSM(1.65547,false,0.000690983),
        		new ScoredPSM(1.65441,false,0.000693639),
        		new ScoredPSM(1.64685,false,0.000712813),
        		new ScoredPSM(1.64161,false,0.000726429),
        		new ScoredPSM(1.63897,false,0.000733385),
        		new ScoredPSM(1.63774,false,0.00073664),
        		new ScoredPSM(1.63145,false,0.000753572),
        		new ScoredPSM(1.63132,false,0.000753932),
        		new ScoredPSM(1.63017,false,0.000757055),
        		new ScoredPSM(1.62982,false,0.00075801),
        		new ScoredPSM(1.62483,false,0.000771789),
        		new ScoredPSM(1.62321,false,0.000776317),
        		new ScoredPSM(1.61668,false,0.000794818),
        		new ScoredPSM(1.6119,false,0.000808653),
        		new ScoredPSM(1.60617,false,0.000825559),
        		new ScoredPSM(1.60381,false,0.000832621),
        		new ScoredPSM(1.60347,false,0.000833654),
        		new ScoredPSM(1.59213,false,0.000868467),
        		new ScoredPSM(1.57519,false,0.00092323),
        		new ScoredPSM(1.57095,false,0.000937437),
        		new ScoredPSM(1.56516,false,0.000957248),
        		new ScoredPSM(1.56305,false,0.000964563),
        		new ScoredPSM(1.55356,false,0.000998143),
        		new ScoredPSM(1.54779,false,0.00101916),
        		new ScoredPSM(1.54155,false,0.00104237),
        		new ScoredPSM(1.53298,false,0.0010751),
        		new ScoredPSM(1.51998,false,0.00112675),
        		new ScoredPSM(1.51864,false,0.0011322),
        		new ScoredPSM(1.50766,false,0.00117795),
        		new ScoredPSM(1.50623,false,0.00118403),
        		new ScoredPSM(1.5044,false,0.00119188),
        		new ScoredPSM(1.49378,false,0.00123843),
        		new ScoredPSM(1.49259,false,0.00124372),
        		new ScoredPSM(1.49117,false,0.00125014),
        		new ScoredPSM(1.48648,false,0.00127147),
        		new ScoredPSM(1.47634,false,0.00131881),
        		new ScoredPSM(1.45942,false,0.00140181),
        		new ScoredPSM(1.45586,false,0.00141992),
        		new ScoredPSM(1.43592,false,0.00152578),
        		new ScoredPSM(1.42041,false,0.00161356),
        		new ScoredPSM(1.41367,false,0.00165323),
        		new ScoredPSM(1.40798,false,0.0016875),
        		new ScoredPSM(1.39016,false,0.00179946),
        		new ScoredPSM(1.38215,false,0.0018522),
        		new ScoredPSM(1.37509,false,0.00189989),
        		new ScoredPSM(1.37253,false,0.00191755),
        		new ScoredPSM(1.37148,false,0.0019248),
        		new ScoredPSM(1.35317,false,0.00205606),
        		new ScoredPSM(1.35137,false,0.00206943),
        		new ScoredPSM(1.34005,false,0.00215561),
        		new ScoredPSM(1.33929,false,0.00216154),
        		new ScoredPSM(1.33292,false,0.00221169),
        		new ScoredPSM(1.33256,false,0.0022146),
        		new ScoredPSM(1.3232,false,0.00229056),
        		new ScoredPSM(1.31083,false,0.00239489),
        		new ScoredPSM(1.30118,false,0.00247957),
        		new ScoredPSM(1.29684,false,0.00251867),
        		new ScoredPSM(1.29631,false,0.00252345),
        		new ScoredPSM(1.29574,false,0.00252867),
        		new ScoredPSM(1.29437,false,0.00254116),
        		new ScoredPSM(1.29067,false,0.00257527),
        		new ScoredPSM(1.27065,false,0.0027677),
        		new ScoredPSM(1.26649,false,0.0028095),
        		new ScoredPSM(1.26433,false,0.00283148),
        		new ScoredPSM(1.25778,false,0.00289894),
        		new ScoredPSM(1.25662,false,0.00291115),
        		new ScoredPSM(1.25384,false,0.00294044),
        		new ScoredPSM(1.24608,false,0.0030237),
        		new ScoredPSM(1.23637,false,0.00313117),
        		new ScoredPSM(1.22884,false,0.00321726),
        		new ScoredPSM(1.22883,false,0.00321737),
        		new ScoredPSM(1.22819,false,0.00322476),
        		new ScoredPSM(1.22244,false,0.00329216),
        		new ScoredPSM(1.20955,false,0.00344838),
        		new ScoredPSM(1.20899,false,0.00345532),
        		new ScoredPSM(1.20708,false,0.0034792),
        		new ScoredPSM(1.20633,false,0.0034886),
        		new ScoredPSM(1.20391,false,0.00351906),
        		new ScoredPSM(1.20061,false,0.00356116),
        		new ScoredPSM(1.18027,false,0.0038313),
        		new ScoredPSM(1.17765,false,0.00386766),
        		new ScoredPSM(1.16949,false,0.00398278),
        		new ScoredPSM(1.16622,false,0.00402983),
        		new ScoredPSM(1.16042,false,0.0041148),
        		new ScoredPSM(1.13221,false,0.00455372),
        		new ScoredPSM(1.10492,false,0.00502247),
        		new ScoredPSM(1.08275,false,0.00543849),
        		new ScoredPSM(1.06514,false,0.00579325),
        		new ScoredPSM(1.05872,false,0.00592811),
        		new ScoredPSM(1.04363,false,0.00625775),
        		new ScoredPSM(1.03969,false,0.00634672),
        		new ScoredPSM(1.00135,false,0.00728094),
        		new ScoredPSM(0.997911,false,0.00737123),
        		new ScoredPSM(0.973455,false,0.00804535),
        		new ScoredPSM(0.965624,false,0.00827385),
        		new ScoredPSM(0.963072,false,0.00834971),
        		new ScoredPSM(0.962096,false,0.0083789),
        		new ScoredPSM(0.955956,false,0.00856485),
        		new ScoredPSM(0.925023,false,0.00956566),
        		new ScoredPSM(0.91427,false,0.00994),
        		new ScoredPSM(0.911362,false,0.0100437),
        		new ScoredPSM(0.909272,false,0.0101189),
        		new ScoredPSM(0.899628,false,0.0104731),
        		new ScoredPSM(0.896952,false,0.0105735),
        		new ScoredPSM(0.866855,false,0.0117709),
        		new ScoredPSM(0.841963,false,0.0128615),
        		new ScoredPSM(0.819314,false,0.0139401),
        		new ScoredPSM(0.810633,false,0.0143767),
        		new ScoredPSM(0.809093,false,0.0144556),
        		new ScoredPSM(0.790665,false,0.0154328),
        		new ScoredPSM(0.772815,false,0.0164413),
        		new ScoredPSM(0.771378,false,0.0165252),
        		new ScoredPSM(0.746335,false,0.0180577),
        		new ScoredPSM(0.740423,false,0.0184394),
        		new ScoredPSM(0.729132,false,0.0191905),
        		new ScoredPSM(0.717274,false,0.0200115),
        		new ScoredPSM(0.716214,false,0.0200866),
        		new ScoredPSM(0.709595,false,0.0205614),
        		new ScoredPSM(0.680541,false,0.0227796),
        		new ScoredPSM(0.675983,false,0.0231481),
        		new ScoredPSM(0.6755,false,0.0231875),
        		new ScoredPSM(0.67386,false,0.0233218),
        		new ScoredPSM(0.671591,false,0.0235088),
        		new ScoredPSM(0.671223,false,0.0235392),
        		new ScoredPSM(0.6674,false,0.0238581),
        		new ScoredPSM(0.656947,false,0.0247515),
        		new ScoredPSM(0.652382,false,0.0251517),
        		new ScoredPSM(0.635214,false,0.0267147),
        		new ScoredPSM(0.58734,false,0.0315887),
        		new ScoredPSM(0.580023,false,0.0324059),
        		new ScoredPSM(0.569764,false,0.0335865),
        		new ScoredPSM(0.562388,false,0.0344608),
        		new ScoredPSM(0.551708,false,0.0357658),
        		new ScoredPSM(0.544147,false,0.0367184),
        		new ScoredPSM(0.524689,false,0.0392841),
        		new ScoredPSM(0.510458,false,0.0412689),
        		new ScoredPSM(0.49741,false,0.0431732),
        		new ScoredPSM(0.436291,false,0.0532699),
        		new ScoredPSM(0.426872,false,0.0550132),
        		new ScoredPSM(0.408509,false,0.058569),
        		new ScoredPSM(0.402302,false,0.0598194),
        		new ScoredPSM(0.402121,false,0.0598564),
        		new ScoredPSM(0.397557,false,0.0607924),
        		new ScoredPSM(0.394076,false,0.0615156),
        		new ScoredPSM(0.371721,false,0.0663569),
        		new ScoredPSM(0.36993,false,0.0667599),
        		new ScoredPSM(0.369379,false,0.0668844),
        		new ScoredPSM(0.360203,false,0.0689889),
        		new ScoredPSM(0.343103,false,0.0730776),
        		new ScoredPSM(0.340181,false,0.0737987),
        		new ScoredPSM(0.333293,false,0.0755243),
        		new ScoredPSM(0.317149,false,0.0797176),
        		new ScoredPSM(0.312647,false,0.0809248),
        		new ScoredPSM(0.299204,false,0.0846315),
        		new ScoredPSM(0.272445,false,0.0924804),
        		new ScoredPSM(0.263029,false,0.0953976),
        		new ScoredPSM(0.244715,false,0.101314),
        		new ScoredPSM(0.242165,false,0.102164),
        		new ScoredPSM(0.231351,false,0.10584),
        		new ScoredPSM(0.228089,false,0.106972),
        		new ScoredPSM(0.223528,false,0.108573),
        		new ScoredPSM(0.222164,false,0.109057),
        		new ScoredPSM(0.221535,false,0.10928),
        		new ScoredPSM(0.201642,false,0.116562),
        		new ScoredPSM(0.199658,false,0.117312),
        		new ScoredPSM(0.183289,false,0.123664),
        		new ScoredPSM(0.172009,false,0.128217),
        		new ScoredPSM(0.166643,false,0.130434),
        		new ScoredPSM(0.165087,false,0.131083),
        		new ScoredPSM(0.158731,false,0.133765),
        		new ScoredPSM(0.140338,false,0.1418),
        		new ScoredPSM(0.13021,false,0.146402),
        		new ScoredPSM(0.123227,false,0.14965),
        		new ScoredPSM(0.115536,false,0.1533),
        		new ScoredPSM(0.110548,false,0.155707),
        		new ScoredPSM(0.0980469,false,0.161883),
        		new ScoredPSM(0.0861998,false,0.167928),
        		new ScoredPSM(0.0605894,false,0.181646),
        		new ScoredPSM(0.0596065,false,0.182191),
        		new ScoredPSM(0.0573393,false,0.183452),
        		new ScoredPSM(0.0572595,false,0.183496),
        		new ScoredPSM(0.0133425,false,0.209378),
        		new ScoredPSM(0.012165,false,0.210111),
        		new ScoredPSM(0,false,0.217798),
        		new ScoredPSM(-0.0046873,false,0.220818),
        		new ScoredPSM(-0.0180205,false,0.22959),
        		new ScoredPSM(-0.0253165,false,0.234504),
        		new ScoredPSM(-0.0302421,false,0.237866),
        		new ScoredPSM(-0.0307954,false,0.238246),
        		new ScoredPSM(-0.0396793,false,0.244412),
        		new ScoredPSM(-0.0458604,false,0.248773),
        		new ScoredPSM(-0.0499802,false,0.251713),
        		new ScoredPSM(-0.0533575,false,0.254142),
        		new ScoredPSM(-0.0550047,false,0.255333),
        		new ScoredPSM(-0.0561518,false,0.256165),
        		new ScoredPSM(-0.0632197,false,0.261335),
        		new ScoredPSM(-0.0637344,false,0.261715),
        		new ScoredPSM(-0.0700351,false,0.266394),
        		new ScoredPSM(-0.076815,false,0.271498),
        		new ScoredPSM(-0.091274,false,0.28262),
        		new ScoredPSM(-0.0948495,false,0.28542),
        		new ScoredPSM(-0.0976583,false,0.287633),
        		new ScoredPSM(-0.10312,false,0.291972),
        		new ScoredPSM(-0.13047,false,0.31439),
        		new ScoredPSM(-0.138643,false,0.321311),
        		new ScoredPSM(-0.145046,false,0.326805),
        		new ScoredPSM(-0.148257,false,0.329583),
        		new ScoredPSM(-0.151166,false,0.332114),
        		new ScoredPSM(-0.156173,false,0.336499),
        		new ScoredPSM(-0.160037,false,0.339909),
        		new ScoredPSM(-0.161514,false,0.341218),
        		new ScoredPSM(-0.165628,false,0.344884),
        		new ScoredPSM(-0.169441,false,0.348303),
        		new ScoredPSM(-0.173654,false,0.352106),
        		new ScoredPSM(-0.177321,false,0.355439),
        		new ScoredPSM(-0.184887,false,0.362376),
        		new ScoredPSM(-0.186895,false,0.364232),
        		new ScoredPSM(-0.201377,false,0.377788),
        		new ScoredPSM(-0.20634,false,0.382504),
        		new ScoredPSM(-0.214786,false,0.39061),
        		new ScoredPSM(-0.221961,false,0.397576),
        		new ScoredPSM(-0.224055,false,0.399622),
        		new ScoredPSM(-0.224436,false,0.399996),
        		new ScoredPSM(-0.227638,false,0.403138),
        		new ScoredPSM(-0.232844,false,0.408278),
        		new ScoredPSM(-0.234723,false,0.410143),
        		new ScoredPSM(-0.235605,false,0.41102),
        		new ScoredPSM(-0.239811,false,0.415216),
        		new ScoredPSM(-0.244824,false,0.420248),
        		new ScoredPSM(-0.248433,false,0.42389),
        		new ScoredPSM(-0.263855,false,0.43965),
        		new ScoredPSM(-0.265563,false,0.441415),
        		new ScoredPSM(-0.267551,false,0.443473),
        		new ScoredPSM(-0.281155,false,0.457685),
        		new ScoredPSM(-0.28429,false,0.460992),
        		new ScoredPSM(-0.288388,false,0.465332),
        		new ScoredPSM(-0.294379,false,0.471711),
        		new ScoredPSM(-0.299525,false,0.477222),
        		new ScoredPSM(-0.311267,false,0.489903),
        		new ScoredPSM(-0.315592,false,0.494609),
        		new ScoredPSM(-0.31996,false,0.499381),
        		new ScoredPSM(-0.321228,false,0.50077),
        		new ScoredPSM(-0.326173,false,0.506202),
        		new ScoredPSM(-0.333694,false,0.514505),
        		new ScoredPSM(-0.361896,false,0.546052),
        		new ScoredPSM(-0.382869,false,0.569863),
        		new ScoredPSM(-0.391828,false,0.580107),
        		new ScoredPSM(-0.39241,false,0.580774),
        		new ScoredPSM(-0.398668,false,0.587953),
        		new ScoredPSM(-0.40463,false,0.594806),
        		new ScoredPSM(-0.427167,false,0.620814),
        		new ScoredPSM(-0.429023,false,0.622962),
        		new ScoredPSM(-0.430755,false,0.624967),
        		new ScoredPSM(-0.432145,false,0.626576),
        		new ScoredPSM(-0.438975,false,0.634487),
        		new ScoredPSM(-0.440289,false,0.636011),
        		new ScoredPSM(-0.445818,false,0.642421),
        		new ScoredPSM(-0.451918,false,0.649495),
        		new ScoredPSM(-0.452177,false,0.649796),
        		new ScoredPSM(-0.46574,false,0.665529),
        		new ScoredPSM(-0.46833,false,0.668532),
        		new ScoredPSM(-0.469957,false,0.67042),
        		new ScoredPSM(-0.47224,false,0.673066),
        		new ScoredPSM(-0.48736,false,0.690577),
        		new ScoredPSM(-0.487859,false,0.691155),
        		new ScoredPSM(-0.488664,false,0.692085),
        		new ScoredPSM(-0.493813,false,0.698037),
        		new ScoredPSM(-0.50117,false,0.706529),
        		new ScoredPSM(-0.507427,false,0.713737),
        		new ScoredPSM(-0.516899,false,0.724624),
        		new ScoredPSM(-0.529154,false,0.738649),
        		new ScoredPSM(-0.530528,false,0.740217),
        		new ScoredPSM(-0.533409,false,0.743502),
        		new ScoredPSM(-0.54468,false,0.756304),
        		new ScoredPSM(-0.545861,false,0.75764),
        		new ScoredPSM(-0.547115,false,0.759059),
        		new ScoredPSM(-0.556681,false,0.769845),
        		new ScoredPSM(-0.561338,false,0.775071),
        		new ScoredPSM(-0.565599,false,0.77984),
        		new ScoredPSM(-0.569323,false,0.783995),
        		new ScoredPSM(-0.569338,false,0.784012),
        		new ScoredPSM(-0.573424,false,0.788558),
        		new ScoredPSM(-0.575156,false,0.790481),
        		new ScoredPSM(-0.578,false,0.793634),
        		new ScoredPSM(-0.583504,false,0.799713),
        		new ScoredPSM(-0.594121,false,0.81136),
        		new ScoredPSM(-0.60267,false,0.820658),
        		new ScoredPSM(-0.616106,false,0.835114),
        		new ScoredPSM(-0.616735,false,0.835787),
        		new ScoredPSM(-0.618315,false,0.837473),
        		new ScoredPSM(-0.620208,false,0.839489),
        		new ScoredPSM(-0.620907,false,0.840233),
        		new ScoredPSM(-0.623424,false,0.842907),
        		new ScoredPSM(-0.626813,false,0.846495),
        		new ScoredPSM(-0.638211,false,0.858465),
        		new ScoredPSM(-0.640124,false,0.860459),
        		new ScoredPSM(-0.643945,false,0.864429),
        		new ScoredPSM(-0.645823,false,0.866373),
        		new ScoredPSM(-0.647555,false,0.868162),
        		new ScoredPSM(-0.647998,false,0.86862),
        		new ScoredPSM(-0.649152,false,0.86981),
        		new ScoredPSM(-0.649287,false,0.869948),
        		new ScoredPSM(-0.656105,false,0.876942),
        		new ScoredPSM(-0.660541,false,0.881461),
        		new ScoredPSM(-0.663307,false,0.884266),
        		new ScoredPSM(-0.671868,false,0.892881),
        		new ScoredPSM(-0.673519,false,0.894531),
        		new ScoredPSM(-0.680951,false,0.901911),
        		new ScoredPSM(-0.682306,false,0.903249),
        		new ScoredPSM(-0.683889,false,0.904808),
        		new ScoredPSM(-0.685475,false,0.906365),
        		new ScoredPSM(-0.689618,false,0.91042),
        		new ScoredPSM(-0.693147,false,0.913852),
        		new ScoredPSM(-0.6933,false,0.914001),
        		new ScoredPSM(-0.693846,false,0.91453),
        		new ScoredPSM(-0.695567,false,0.916195),
        		new ScoredPSM(-0.695664,false,0.916289),
        		new ScoredPSM(-0.698375,false,0.918903),
        		new ScoredPSM(-0.704181,false,0.924464),
        		new ScoredPSM(-0.714648,false,0.934351),
        		new ScoredPSM(-0.716975,false,0.936525),
        		new ScoredPSM(-0.725295,false,0.944223),
        		new ScoredPSM(-0.728711,false,0.947349),
        		new ScoredPSM(-0.729646,false,0.948201),
        		new ScoredPSM(-0.739268,false,0.956878),
        		new ScoredPSM(-0.739528,false,0.95711),
        		new ScoredPSM(-0.742638,false,0.959878),
        		new ScoredPSM(-0.754015,false,0.969848),
        		new ScoredPSM(-0.756791,false,0.972244),
        		new ScoredPSM(-0.757013,false,0.972435),
        		new ScoredPSM(-0.758202,false,0.973455),
        		new ScoredPSM(-0.758382,false,0.97361),
        		new ScoredPSM(-0.765466,false,0.979634),
        		new ScoredPSM(-0.767103,false,0.981012),
        		new ScoredPSM(-0.770943,false,0.984224),
        		new ScoredPSM(-0.786016,false,0.996545),
        		new ScoredPSM(-0.793041,false,1),
        		new ScoredPSM(-0.798784,false,1),
        		new ScoredPSM(-0.806354,false,1),
        		new ScoredPSM(-0.813329,false,1),
        		new ScoredPSM(-0.826842,false,1),
        		new ScoredPSM(-0.828257,false,1),
        		new ScoredPSM(-0.838823,false,1),
        		new ScoredPSM(-0.845716,false,1),
        		new ScoredPSM(-0.846017,false,1),
        		new ScoredPSM(-0.850115,false,1),
        		new ScoredPSM(-0.851195,false,1),
        		new ScoredPSM(-0.854291,false,1),
        		new ScoredPSM(-0.856292,false,1),
        		new ScoredPSM(-0.857613,false,1),
        		new ScoredPSM(-0.860589,false,1),
        		new ScoredPSM(-0.888563,false,1),
        		new ScoredPSM(-0.889984,false,1),
        		new ScoredPSM(-0.893799,false,1),
        		new ScoredPSM(-0.902901,false,1),
        		new ScoredPSM(-0.903742,false,1),
        		new ScoredPSM(-0.90652,false,1),
        		new ScoredPSM(-0.910665,false,1),
        		new ScoredPSM(-0.914346,false,1),
        		new ScoredPSM(-0.91972,false,1),
        		new ScoredPSM(-0.934284,false,1),
        		new ScoredPSM(-0.942464,false,1),
        		new ScoredPSM(-0.951685,false,1),
        		new ScoredPSM(-0.952569,false,1),
        		new ScoredPSM(-0.957165,false,1),
        		new ScoredPSM(-0.957528,false,1),
        		new ScoredPSM(-0.959374,false,1),
        		new ScoredPSM(-0.95954,false,1),
        		new ScoredPSM(-0.961758,false,1),
        		new ScoredPSM(-0.964401,false,1),
        		new ScoredPSM(-0.96586,false,1),
        		new ScoredPSM(-0.969386,false,1),
        		new ScoredPSM(-0.970071,false,1),
        		new ScoredPSM(-0.970503,false,1),
        		new ScoredPSM(-0.984052,false,1),
        		new ScoredPSM(-0.988015,false,1),
        		new ScoredPSM(-1.00329,false,1),
        		new ScoredPSM(-1.0034,false,1),
        		new ScoredPSM(-1.00821,false,1),
        		new ScoredPSM(-1.01139,false,1),
        		new ScoredPSM(-1.0159,false,1),
        		new ScoredPSM(-1.01986,false,1),
        		new ScoredPSM(-1.02988,false,1),
        		new ScoredPSM(-1.03241,false,1),
        		new ScoredPSM(-1.0349,false,1),
        		new ScoredPSM(-1.04247,false,1),
        		new ScoredPSM(-1.04474,false,1),
        		new ScoredPSM(-1.04885,false,1),
        		new ScoredPSM(-1.05047,false,1),
        		new ScoredPSM(-1.06119,false,1),
        		new ScoredPSM(-1.06645,false,1),
        		new ScoredPSM(-1.07258,false,1),
        		new ScoredPSM(-1.07618,false,1),
        		new ScoredPSM(-1.08117,false,1),
        		new ScoredPSM(-1.09204,false,1),
        		new ScoredPSM(-1.09653,false,1),
        		new ScoredPSM(-1.09967,false,1),
        		new ScoredPSM(-1.11173,false,1),
        		new ScoredPSM(-1.12136,false,1),
        		new ScoredPSM(-1.12805,false,1),
        		new ScoredPSM(-1.13139,false,1),
        		new ScoredPSM(-1.13755,false,1),
        		new ScoredPSM(-1.13805,false,1),
        		new ScoredPSM(-1.14203,false,1),
        		new ScoredPSM(-1.14585,false,1),
        		new ScoredPSM(-1.15669,false,1),
        		new ScoredPSM(-1.15845,false,1),
        		new ScoredPSM(-1.16052,false,1),
        		new ScoredPSM(-1.16278,false,1),
        		new ScoredPSM(-1.16395,false,1),
        		new ScoredPSM(-1.17644,false,1),
        		new ScoredPSM(-1.17884,false,1),
        		new ScoredPSM(-1.17887,false,1),
        		new ScoredPSM(-1.18208,false,1),
        		new ScoredPSM(-1.18887,false,1),
        		new ScoredPSM(-1.19673,false,1),
        		new ScoredPSM(-1.21896,false,1),
        		new ScoredPSM(-1.21949,false,1),
        		new ScoredPSM(-1.22643,false,1),
        		new ScoredPSM(-1.23934,false,1),
        		new ScoredPSM(-1.23971,false,1),
        		new ScoredPSM(-1.24025,false,1),
        		new ScoredPSM(-1.24557,false,1),
        		new ScoredPSM(-1.24886,false,1),
        		new ScoredPSM(-1.25087,false,1),
        		new ScoredPSM(-1.25157,false,1),
        		new ScoredPSM(-1.25997,false,1),
        		new ScoredPSM(-1.27378,false,1),
        		new ScoredPSM(-1.27413,false,1),
        		new ScoredPSM(-1.28652,false,1),
        		new ScoredPSM(-1.29698,false,1),
        		new ScoredPSM(-1.29738,false,1),
        		new ScoredPSM(-1.30869,false,1),
        		new ScoredPSM(-1.31484,false,1),
        		new ScoredPSM(-1.33033,false,1),
        		new ScoredPSM(-1.3331,false,1),
        		new ScoredPSM(-1.3336,false,1),
        		new ScoredPSM(-1.33958,false,1),
        		new ScoredPSM(-1.3518,false,1),
        		new ScoredPSM(-1.35801,false,1),
        		new ScoredPSM(-1.36839,false,1),
        		new ScoredPSM(-1.37257,false,1),
        		new ScoredPSM(-1.37265,false,1),
        		new ScoredPSM(-1.37314,false,1),
        		new ScoredPSM(-1.38298,false,1),
        		new ScoredPSM(-1.38989,false,1),
        		new ScoredPSM(-1.40335,false,1),
        		new ScoredPSM(-1.40687,false,1),
        		new ScoredPSM(-1.40847,false,1),
        		new ScoredPSM(-1.41756,false,1),
        		new ScoredPSM(-1.41916,false,1),
        		new ScoredPSM(-1.43001,false,1),
        		new ScoredPSM(-1.43566,false,1),
        		new ScoredPSM(-1.44808,false,1),
        		new ScoredPSM(-1.44838,false,1),
        		new ScoredPSM(-1.44998,false,1),
        		new ScoredPSM(-1.45796,false,1),
        		new ScoredPSM(-1.45939,false,1),
        		new ScoredPSM(-1.47484,false,1),
        		new ScoredPSM(-1.47576,false,1),
        		new ScoredPSM(-1.51097,false,1),
        		new ScoredPSM(-1.52877,false,1),
        		new ScoredPSM(-1.54934,false,1),
        		new ScoredPSM(-1.5501,false,1),
        		new ScoredPSM(-1.58001,false,1),
        		new ScoredPSM(-1.58119,false,1),
        		new ScoredPSM(-1.58238,false,1),
        		new ScoredPSM(-1.60103,false,1),
        		new ScoredPSM(-1.60344,false,1),
        		new ScoredPSM(-1.61183,false,1),
        		new ScoredPSM(-1.62867,false,1),
        		new ScoredPSM(-1.63043,false,1),
        		new ScoredPSM(-1.6314,false,1),
        		new ScoredPSM(-1.63757,false,1),
        		new ScoredPSM(-1.65395,false,1),
        		new ScoredPSM(-1.65552,false,1),
        		new ScoredPSM(-1.66277,false,1),
        		new ScoredPSM(-1.67248,false,1),
        		new ScoredPSM(-1.6791,false,1),
        		new ScoredPSM(-1.68213,false,1),
        		new ScoredPSM(-1.71678,false,1),
        		new ScoredPSM(-1.7212,false,1),
        		new ScoredPSM(-1.76874,false,1),
        		new ScoredPSM(-1.7693,false,1),
        		new ScoredPSM(-1.78047,false,1),
        		new ScoredPSM(-1.80638,false,1),
        		new ScoredPSM(-1.82748,false,1),
        		new ScoredPSM(-1.84294,false,1),
        		new ScoredPSM(-1.95346,false,1),
        		new ScoredPSM(-1.95919,false,1),
        		new ScoredPSM(-1.95973,false,1),
        		new ScoredPSM(-2.02352,false,1),
        		new ScoredPSM(-2.12829,false,1),
        		new ScoredPSM(-2.16962,false,1),
        		new ScoredPSM(-2.17086,false,1),
        		new ScoredPSM(-2.20847,false,1),
        		new ScoredPSM(-2.28883,false,1),
        		new ScoredPSM(-2.34421,false,1),
        		new ScoredPSM(-2.46027,false,1),
        		new ScoredPSM(-2.47455,false,1),
        		new ScoredPSM(-2.50975,false,1),
        		new ScoredPSM(-2.54552,false,1),
        		new ScoredPSM(-2.56062,false,1),
        		new ScoredPSM(-2.77573,false,1),
        		new ScoredPSM(0.348682,true,0.0717195),
        		new ScoredPSM(0.212267,true,0.11262),
        		new ScoredPSM(0.191744,true,0.120346),
        		new ScoredPSM(0.0784242,true,0.171998),
        		new ScoredPSM(0.0710491,true,0.175934),
        		new ScoredPSM(0.0356396,true,0.195891),
        		new ScoredPSM(-0.00174875,true,0.218921),
        		new ScoredPSM(-0.0142079,true,0.227055),
        		new ScoredPSM(-0.01522,true,0.227726),
        		new ScoredPSM(-0.0188774,true,0.230163),
        		new ScoredPSM(-0.0456091,true,0.248595),
        		new ScoredPSM(-0.0460536,true,0.248911),
        		new ScoredPSM(-0.0526602,true,0.253639),
        		new ScoredPSM(-0.0625233,true,0.260822),
        		new ScoredPSM(-0.0792697,true,0.273363),
        		new ScoredPSM(-0.0794661,true,0.273513),
        		new ScoredPSM(-0.104874,true,0.293375),
        		new ScoredPSM(-0.11157,true,0.298776),
        		new ScoredPSM(-0.114042,true,0.300787),
        		new ScoredPSM(-0.120253,true,0.305882),
        		new ScoredPSM(-0.126191,true,0.310807),
        		new ScoredPSM(-0.132208,true,0.315853),
        		new ScoredPSM(-0.132528,true,0.316123),
        		new ScoredPSM(-0.135393,true,0.318547),
        		new ScoredPSM(-0.158231,true,0.338312),
        		new ScoredPSM(-0.1586,true,0.338639),
        		new ScoredPSM(-0.162496,true,0.342091),
        		new ScoredPSM(-0.170879,true,0.349598),
        		new ScoredPSM(-0.17317,true,0.351668),
        		new ScoredPSM(-0.179797,true,0.3577),
        		new ScoredPSM(-0.1873,true,0.364607),
        		new ScoredPSM(-0.190304,true,0.367395),
        		new ScoredPSM(-0.197568,true,0.374193),
        		new ScoredPSM(-0.210042,true,0.386044),
        		new ScoredPSM(-0.215779,true,0.39157),
        		new ScoredPSM(-0.229269,true,0.404744),
        		new ScoredPSM(-0.230169,true,0.405633),
        		new ScoredPSM(-0.232332,true,0.407771),
        		new ScoredPSM(-0.236197,true,0.411608),
        		new ScoredPSM(-0.239933,true,0.415337),
        		new ScoredPSM(-0.245457,true,0.420885),
        		new ScoredPSM(-0.264966,true,0.440798),
        		new ScoredPSM(-0.279637,true,0.456088),
        		new ScoredPSM(-0.28324,true,0.459883),
        		new ScoredPSM(-0.291465,true,0.468603),
        		new ScoredPSM(-0.292724,true,0.469944),
        		new ScoredPSM(-0.299338,true,0.477021),
        		new ScoredPSM(-0.309261,true,0.487726),
        		new ScoredPSM(-0.310321,true,0.488875),
        		new ScoredPSM(-0.316424,true,0.495517),
        		new ScoredPSM(-0.319256,true,0.498611),
        		new ScoredPSM(-0.323105,true,0.502829),
        		new ScoredPSM(-0.331268,true,0.511821),
        		new ScoredPSM(-0.333756,true,0.514574),
        		new ScoredPSM(-0.336855,true,0.518009),
        		new ScoredPSM(-0.350635,true,0.533383),
        		new ScoredPSM(-0.35426,true,0.537451),
        		new ScoredPSM(-0.355752,true,0.539129),
        		new ScoredPSM(-0.357451,true,0.541041),
        		new ScoredPSM(-0.367451,true,0.552333),
        		new ScoredPSM(-0.371076,true,0.556442),
        		new ScoredPSM(-0.372781,true,0.558378),
        		new ScoredPSM(-0.373184,true,0.558835),
        		new ScoredPSM(-0.375023,true,0.560926),
        		new ScoredPSM(-0.377975,true,0.564285),
        		new ScoredPSM(-0.388555,true,0.57636),
        		new ScoredPSM(-0.391601,true,0.579847),
        		new ScoredPSM(-0.397489,true,0.586599),
        		new ScoredPSM(-0.397998,true,0.587183),
        		new ScoredPSM(-0.401907,true,0.591674),
        		new ScoredPSM(-0.402933,true,0.592853),
        		new ScoredPSM(-0.40358,true,0.593598),
        		new ScoredPSM(-0.405792,true,0.596143),
        		new ScoredPSM(-0.406624,true,0.597101),
        		new ScoredPSM(-0.406879,true,0.597394),
        		new ScoredPSM(-0.411599,true,0.602833),
        		new ScoredPSM(-0.412031,true,0.603331),
        		new ScoredPSM(-0.41908,true,0.611466),
        		new ScoredPSM(-0.425926,true,0.619379),
        		new ScoredPSM(-0.426079,true,0.619556),
        		new ScoredPSM(-0.427915,true,0.62168),
        		new ScoredPSM(-0.428362,true,0.622197),
        		new ScoredPSM(-0.428818,true,0.622725),
        		new ScoredPSM(-0.430704,true,0.624908),
        		new ScoredPSM(-0.431933,true,0.62633),
        		new ScoredPSM(-0.441312,true,0.637196),
        		new ScoredPSM(-0.445601,true,0.642169),
        		new ScoredPSM(-0.453002,true,0.650752),
        		new ScoredPSM(-0.46181,true,0.66097),
        		new ScoredPSM(-0.46343,true,0.66285),
        		new ScoredPSM(-0.464314,true,0.663876),
        		new ScoredPSM(-0.466765,true,0.666718),
        		new ScoredPSM(-0.469972,true,0.670437),
        		new ScoredPSM(-0.471915,true,0.672689),
        		new ScoredPSM(-0.472224,true,0.673047),
        		new ScoredPSM(-0.474136,true,0.675264),
        		new ScoredPSM(-0.474148,true,0.675277),
        		new ScoredPSM(-0.474481,true,0.675664),
        		new ScoredPSM(-0.474512,true,0.6757),
        		new ScoredPSM(-0.478186,true,0.679957),
        		new ScoredPSM(-0.484471,true,0.687235),
        		new ScoredPSM(-0.484984,true,0.687829),
        		new ScoredPSM(-0.485086,true,0.687946),
        		new ScoredPSM(-0.488119,true,0.691456),
        		new ScoredPSM(-0.495728,true,0.700249),
        		new ScoredPSM(-0.508431,true,0.714893),
        		new ScoredPSM(-0.512145,true,0.719165),
        		new ScoredPSM(-0.522619,true,0.731179),
        		new ScoredPSM(-0.526232,true,0.735312),
        		new ScoredPSM(-0.526599,true,0.735731),
        		new ScoredPSM(-0.526641,true,0.735779),
        		new ScoredPSM(-0.529559,true,0.739111),
        		new ScoredPSM(-0.531613,true,0.741455),
        		new ScoredPSM(-0.532888,true,0.742908),
        		new ScoredPSM(-0.538605,true,0.749412),
        		new ScoredPSM(-0.538943,true,0.749796),
        		new ScoredPSM(-0.542864,true,0.754245),
        		new ScoredPSM(-0.544164,true,0.755719),
        		new ScoredPSM(-0.545836,true,0.757612),
        		new ScoredPSM(-0.546665,true,0.75855),
        		new ScoredPSM(-0.551007,true,0.763455),
        		new ScoredPSM(-0.551026,true,0.763476),
        		new ScoredPSM(-0.552313,true,0.764927),
        		new ScoredPSM(-0.55389,true,0.766704),
        		new ScoredPSM(-0.555505,true,0.768522),
        		new ScoredPSM(-0.556723,true,0.769892),
        		new ScoredPSM(-0.559705,true,0.77324),
        		new ScoredPSM(-0.560782,true,0.774448),
        		new ScoredPSM(-0.561351,true,0.775086),
        		new ScoredPSM(-0.563024,true,0.776959),
        		new ScoredPSM(-0.569205,true,0.783864),
        		new ScoredPSM(-0.56964,true,0.784348),
        		new ScoredPSM(-0.571605,true,0.786536),
        		new ScoredPSM(-0.572209,true,0.787207),
        		new ScoredPSM(-0.574624,true,0.789891),
        		new ScoredPSM(-0.578596,true,0.794293),
        		new ScoredPSM(-0.579422,true,0.795207),
        		new ScoredPSM(-0.586623,true,0.803145),
        		new ScoredPSM(-0.58718,true,0.803758),
        		new ScoredPSM(-0.588733,true,0.805463),
        		new ScoredPSM(-0.589813,true,0.806647),
        		new ScoredPSM(-0.590169,true,0.807038),
        		new ScoredPSM(-0.591284,true,0.808258),
        		new ScoredPSM(-0.59338,true,0.810551),
        		new ScoredPSM(-0.594549,true,0.811827),
        		new ScoredPSM(-0.595142,true,0.812474),
        		new ScoredPSM(-0.595516,true,0.812882),
        		new ScoredPSM(-0.596119,true,0.81354),
        		new ScoredPSM(-0.598108,true,0.815705),
        		new ScoredPSM(-0.598153,true,0.815755),
        		new ScoredPSM(-0.598342,true,0.81596),
        		new ScoredPSM(-0.598741,true,0.816394),
        		new ScoredPSM(-0.599696,true,0.817432),
        		new ScoredPSM(-0.600802,true,0.818633),
        		new ScoredPSM(-0.601356,true,0.819233),
        		new ScoredPSM(-0.606493,true,0.82479),
        		new ScoredPSM(-0.606908,true,0.825239),
        		new ScoredPSM(-0.611708,true,0.830404),
        		new ScoredPSM(-0.612069,true,0.830792),
        		new ScoredPSM(-0.612325,true,0.831066),
        		new ScoredPSM(-0.613714,true,0.832555),
        		new ScoredPSM(-0.618717,true,0.837902),
        		new ScoredPSM(-0.620452,true,0.839749),
        		new ScoredPSM(-0.622153,true,0.841557),
        		new ScoredPSM(-0.622727,true,0.842167),
        		new ScoredPSM(-0.624974,true,0.844549),
        		new ScoredPSM(-0.62775,true,0.847484),
        		new ScoredPSM(-0.634178,true,0.854247),
        		new ScoredPSM(-0.634678,true,0.85477),
        		new ScoredPSM(-0.637067,true,0.85727),
        		new ScoredPSM(-0.637801,true,0.858037),
        		new ScoredPSM(-0.639782,true,0.860103),
        		new ScoredPSM(-0.64194,true,0.862348),
        		new ScoredPSM(-0.6432,true,0.863656),
        		new ScoredPSM(-0.645299,true,0.865831),
        		new ScoredPSM(-0.647372,true,0.867973),
        		new ScoredPSM(-0.65355,true,0.874328),
        		new ScoredPSM(-0.653672,true,0.874453),
        		new ScoredPSM(-0.653941,true,0.874729),
        		new ScoredPSM(-0.656103,true,0.87694),
        		new ScoredPSM(-0.657758,true,0.878629),
        		new ScoredPSM(-0.662902,true,0.883855),
        		new ScoredPSM(-0.664636,true,0.885609),
        		new ScoredPSM(-0.66639,true,0.887379),
        		new ScoredPSM(-0.666726,true,0.887718),
        		new ScoredPSM(-0.669487,true,0.890494),
        		new ScoredPSM(-0.669901,true,0.89091),
        		new ScoredPSM(-0.671071,true,0.892083),
        		new ScoredPSM(-0.673603,true,0.894615),
        		new ScoredPSM(-0.679832,true,0.900805),
        		new ScoredPSM(-0.681653,true,0.902604),
        		new ScoredPSM(-0.68298,true,0.903913),
        		new ScoredPSM(-0.685092,true,0.90599),
        		new ScoredPSM(-0.685122,true,0.906019),
        		new ScoredPSM(-0.68661,true,0.907479),
        		new ScoredPSM(-0.689477,true,0.910282),
        		new ScoredPSM(-0.689946,true,0.910739),
        		new ScoredPSM(-0.69128,true,0.912038),
        		new ScoredPSM(-0.69203,true,0.912768),
        		new ScoredPSM(-0.693202,true,0.913905),
        		new ScoredPSM(-0.694084,true,0.914761),
        		new ScoredPSM(-0.698345,true,0.918875),
        		new ScoredPSM(-0.700368,true,0.920817),
        		new ScoredPSM(-0.702625,true,0.922979),
        		new ScoredPSM(-0.705829,true,0.926033),
        		new ScoredPSM(-0.70596,true,0.926157),
        		new ScoredPSM(-0.707626,true,0.927737),
        		new ScoredPSM(-0.707733,true,0.927839),
        		new ScoredPSM(-0.709236,true,0.929261),
        		new ScoredPSM(-0.711283,true,0.931193),
        		new ScoredPSM(-0.711726,true,0.93161),
        		new ScoredPSM(-0.713535,true,0.933309),
        		new ScoredPSM(-0.715365,true,0.935023),
        		new ScoredPSM(-0.715777,true,0.935407),
        		new ScoredPSM(-0.716233,true,0.935833),
        		new ScoredPSM(-0.716286,true,0.935883),
        		new ScoredPSM(-0.717835,true,0.937326),
        		new ScoredPSM(-0.719681,true,0.939042),
        		new ScoredPSM(-0.722514,true,0.941663),
        		new ScoredPSM(-0.726176,true,0.945031),
        		new ScoredPSM(-0.727448,true,0.946196),
        		new ScoredPSM(-0.728823,true,0.947451),
        		new ScoredPSM(-0.731605,true,0.949981),
        		new ScoredPSM(-0.732634,true,0.950913),
        		new ScoredPSM(-0.7328,true,0.951064),
        		new ScoredPSM(-0.734837,true,0.952903),
        		new ScoredPSM(-0.735578,true,0.953571),
        		new ScoredPSM(-0.73685,true,0.954713),
        		new ScoredPSM(-0.737651,true,0.955432),
        		new ScoredPSM(-0.740008,true,0.957539),
        		new ScoredPSM(-0.74391,true,0.961005),
        		new ScoredPSM(-0.745528,true,0.962434),
        		new ScoredPSM(-0.746108,true,0.962945),
        		new ScoredPSM(-0.746905,true,0.963646),
        		new ScoredPSM(-0.746938,true,0.963675),
        		new ScoredPSM(-0.748005,true,0.964612),
        		new ScoredPSM(-0.748403,true,0.96496),
        		new ScoredPSM(-0.748541,true,0.965081),
        		new ScoredPSM(-0.749844,true,0.966221),
        		new ScoredPSM(-0.751658,true,0.967802),
        		new ScoredPSM(-0.751746,true,0.96788),
        		new ScoredPSM(-0.753738,true,0.969608),
        		new ScoredPSM(-0.754823,true,0.970547),
        		new ScoredPSM(-0.755496,true,0.971128),
        		new ScoredPSM(-0.755867,true,0.971448),
        		new ScoredPSM(-0.755927,true,0.9715),
        		new ScoredPSM(-0.760089,true,0.97507),
        		new ScoredPSM(-0.763165,true,0.977688),
        		new ScoredPSM(-0.763886,true,0.978299),
        		new ScoredPSM(-0.767009,true,0.980933),
        		new ScoredPSM(-0.76728,true,0.981161),
        		new ScoredPSM(-0.769355,true,0.9829),
        		new ScoredPSM(-0.774496,true,0.98717),
        		new ScoredPSM(-0.7774,true,0.989559),
        		new ScoredPSM(-0.779051,true,0.990909),
        		new ScoredPSM(-0.779796,true,0.991517),
        		new ScoredPSM(-0.78008,true,0.991748),
        		new ScoredPSM(-0.780701,true,0.992253),
        		new ScoredPSM(-0.78437,true,0.995222),
        		new ScoredPSM(-0.787586,true,0.997801),
        		new ScoredPSM(-0.792754,true,1),
        		new ScoredPSM(-0.793837,true,1),
        		new ScoredPSM(-0.794284,true,1),
        		new ScoredPSM(-0.795855,true,1),
        		new ScoredPSM(-0.803536,true,1),
        		new ScoredPSM(-0.803641,true,1),
        		new ScoredPSM(-0.804637,true,1),
        		new ScoredPSM(-0.804664,true,1),
        		new ScoredPSM(-0.804746,true,1),
        		new ScoredPSM(-0.804891,true,1),
        		new ScoredPSM(-0.807273,true,1),
        		new ScoredPSM(-0.811347,true,1),
        		new ScoredPSM(-0.813913,true,1),
        		new ScoredPSM(-0.816811,true,1),
        		new ScoredPSM(-0.818398,true,1),
        		new ScoredPSM(-0.819347,true,1),
        		new ScoredPSM(-0.820686,true,1),
        		new ScoredPSM(-0.82079,true,1),
        		new ScoredPSM(-0.821046,true,1),
        		new ScoredPSM(-0.825052,true,1),
        		new ScoredPSM(-0.826311,true,1),
        		new ScoredPSM(-0.828057,true,1),
        		new ScoredPSM(-0.828088,true,1),
        		new ScoredPSM(-0.829301,true,1),
        		new ScoredPSM(-0.829835,true,1),
        		new ScoredPSM(-0.831087,true,1),
        		new ScoredPSM(-0.832129,true,1),
        		new ScoredPSM(-0.83632,true,1),
        		new ScoredPSM(-0.838928,true,1),
        		new ScoredPSM(-0.84297,true,1),
        		new ScoredPSM(-0.843515,true,1),
        		new ScoredPSM(-0.844208,true,1),
        		new ScoredPSM(-0.846112,true,1),
        		new ScoredPSM(-0.846983,true,1),
        		new ScoredPSM(-0.85054,true,1),
        		new ScoredPSM(-0.851038,true,1),
        		new ScoredPSM(-0.85148,true,1),
        		new ScoredPSM(-0.851682,true,1),
        		new ScoredPSM(-0.857576,true,1),
        		new ScoredPSM(-0.860725,true,1),
        		new ScoredPSM(-0.862362,true,1),
        		new ScoredPSM(-0.865022,true,1),
        		new ScoredPSM(-0.865735,true,1),
        		new ScoredPSM(-0.866815,true,1),
        		new ScoredPSM(-0.869556,true,1),
        		new ScoredPSM(-0.870781,true,1),
        		new ScoredPSM(-0.870866,true,1),
        		new ScoredPSM(-0.874948,true,1),
        		new ScoredPSM(-0.878114,true,1),
        		new ScoredPSM(-0.878347,true,1),
        		new ScoredPSM(-0.878398,true,1),
        		new ScoredPSM(-0.878967,true,1),
        		new ScoredPSM(-0.881963,true,1),
        		new ScoredPSM(-0.883916,true,1),
        		new ScoredPSM(-0.885093,true,1),
        		new ScoredPSM(-0.885865,true,1),
        		new ScoredPSM(-0.887571,true,1),
        		new ScoredPSM(-0.889481,true,1),
        		new ScoredPSM(-0.894212,true,1),
        		new ScoredPSM(-0.897527,true,1),
        		new ScoredPSM(-0.90384,true,1),
        		new ScoredPSM(-0.9059,true,1),
        		new ScoredPSM(-0.90784,true,1),
        		new ScoredPSM(-0.909327,true,1),
        		new ScoredPSM(-0.910092,true,1),
        		new ScoredPSM(-0.910132,true,1),
        		new ScoredPSM(-0.911958,true,1),
        		new ScoredPSM(-0.912431,true,1),
        		new ScoredPSM(-0.91386,true,1),
        		new ScoredPSM(-0.920251,true,1),
        		new ScoredPSM(-0.920662,true,1),
        		new ScoredPSM(-0.922618,true,1),
        		new ScoredPSM(-0.922722,true,1),
        		new ScoredPSM(-0.923532,true,1),
        		new ScoredPSM(-0.927761,true,1),
        		new ScoredPSM(-0.928098,true,1),
        		new ScoredPSM(-0.928142,true,1),
        		new ScoredPSM(-0.929163,true,1),
        		new ScoredPSM(-0.929598,true,1),
        		new ScoredPSM(-0.935729,true,1),
        		new ScoredPSM(-0.936927,true,1),
        		new ScoredPSM(-0.938826,true,1),
        		new ScoredPSM(-0.940132,true,1),
        		new ScoredPSM(-0.941279,true,1),
        		new ScoredPSM(-0.943557,true,1),
        		new ScoredPSM(-0.943629,true,1),
        		new ScoredPSM(-0.943854,true,1),
        		new ScoredPSM(-0.943901,true,1),
        		new ScoredPSM(-0.946078,true,1),
        		new ScoredPSM(-0.946078,true,1),
        		new ScoredPSM(-0.948748,true,1),
        		new ScoredPSM(-0.951256,true,1),
        		new ScoredPSM(-0.953029,true,1),
        		new ScoredPSM(-0.955137,true,1),
        		new ScoredPSM(-0.956424,true,1),
        		new ScoredPSM(-0.958035,true,1),
        		new ScoredPSM(-0.959005,true,1),
        		new ScoredPSM(-0.962433,true,1),
        		new ScoredPSM(-0.968701,true,1),
        		new ScoredPSM(-0.968708,true,1),
        		new ScoredPSM(-0.969356,true,1),
        		new ScoredPSM(-0.970916,true,1),
        		new ScoredPSM(-0.971886,true,1),
        		new ScoredPSM(-0.979878,true,1),
        		new ScoredPSM(-0.980707,true,1),
        		new ScoredPSM(-0.980925,true,1),
        		new ScoredPSM(-0.984452,true,1),
        		new ScoredPSM(-0.98541,true,1),
        		new ScoredPSM(-0.985627,true,1),
        		new ScoredPSM(-0.989097,true,1),
        		new ScoredPSM(-0.989173,true,1),
        		new ScoredPSM(-0.991036,true,1),
        		new ScoredPSM(-0.995484,true,1),
        		new ScoredPSM(-0.997612,true,1),
        		new ScoredPSM(-0.998719,true,1),
        		new ScoredPSM(-0.998937,true,1),
        		new ScoredPSM(-1,true,1),
        		new ScoredPSM(-1.00293,true,1),
        		new ScoredPSM(-1.00352,true,1),
        		new ScoredPSM(-1.00714,true,1),
        		new ScoredPSM(-1.00776,true,1),
        		new ScoredPSM(-1.00899,true,1),
        		new ScoredPSM(-1.01037,true,1),
        		new ScoredPSM(-1.01334,true,1),
        		new ScoredPSM(-1.01384,true,1),
        		new ScoredPSM(-1.01467,true,1),
        		new ScoredPSM(-1.01669,true,1),
        		new ScoredPSM(-1.02138,true,1),
        		new ScoredPSM(-1.0228,true,1),
        		new ScoredPSM(-1.0237,true,1),
        		new ScoredPSM(-1.0237,true,1),
        		new ScoredPSM(-1.03105,true,1),
        		new ScoredPSM(-1.03266,true,1),
        		new ScoredPSM(-1.03457,true,1),
        		new ScoredPSM(-1.03674,true,1),
        		new ScoredPSM(-1.04015,true,1),
        		new ScoredPSM(-1.04034,true,1),
        		new ScoredPSM(-1.04102,true,1),
        		new ScoredPSM(-1.04113,true,1),
        		new ScoredPSM(-1.04538,true,1),
        		new ScoredPSM(-1.0465,true,1),
        		new ScoredPSM(-1.04789,true,1),
        		new ScoredPSM(-1.05575,true,1),
        		new ScoredPSM(-1.05677,true,1),
        		new ScoredPSM(-1.05782,true,1),
        		new ScoredPSM(-1.06025,true,1),
        		new ScoredPSM(-1.06364,true,1),
        		new ScoredPSM(-1.06512,true,1),
        		new ScoredPSM(-1.06515,true,1),
        		new ScoredPSM(-1.06676,true,1),
        		new ScoredPSM(-1.07191,true,1),
        		new ScoredPSM(-1.07348,true,1),
        		new ScoredPSM(-1.0738,true,1),
        		new ScoredPSM(-1.07521,true,1),
        		new ScoredPSM(-1.0773,true,1),
        		new ScoredPSM(-1.07775,true,1),
        		new ScoredPSM(-1.07842,true,1),
        		new ScoredPSM(-1.07877,true,1),
        		new ScoredPSM(-1.07969,true,1),
        		new ScoredPSM(-1.08089,true,1),
        		new ScoredPSM(-1.08327,true,1),
        		new ScoredPSM(-1.08826,true,1),
        		new ScoredPSM(-1.09158,true,1),
        		new ScoredPSM(-1.09369,true,1),
        		new ScoredPSM(-1.09809,true,1),
        		new ScoredPSM(-1.09941,true,1),
        		new ScoredPSM(-1.09986,true,1),
        		new ScoredPSM(-1.10241,true,1),
        		new ScoredPSM(-1.10315,true,1),
        		new ScoredPSM(-1.10534,true,1),
        		new ScoredPSM(-1.11081,true,1),
        		new ScoredPSM(-1.11218,true,1),
        		new ScoredPSM(-1.11235,true,1),
        		new ScoredPSM(-1.11358,true,1),
        		new ScoredPSM(-1.11391,true,1),
        		new ScoredPSM(-1.11492,true,1),
        		new ScoredPSM(-1.1163,true,1),
        		new ScoredPSM(-1.12015,true,1),
        		new ScoredPSM(-1.12097,true,1),
        		new ScoredPSM(-1.1233,true,1),
        		new ScoredPSM(-1.12347,true,1),
        		new ScoredPSM(-1.12882,true,1),
        		new ScoredPSM(-1.12977,true,1),
        		new ScoredPSM(-1.13052,true,1),
        		new ScoredPSM(-1.13064,true,1),
        		new ScoredPSM(-1.13238,true,1),
        		new ScoredPSM(-1.13328,true,1),
        		new ScoredPSM(-1.1342,true,1),
        		new ScoredPSM(-1.13626,true,1),
        		new ScoredPSM(-1.13662,true,1),
        		new ScoredPSM(-1.13962,true,1),
        		new ScoredPSM(-1.14004,true,1),
        		new ScoredPSM(-1.14149,true,1),
        		new ScoredPSM(-1.14164,true,1),
        		new ScoredPSM(-1.14213,true,1),
        		new ScoredPSM(-1.14749,true,1),
        		new ScoredPSM(-1.14996,true,1),
        		new ScoredPSM(-1.15153,true,1),
        		new ScoredPSM(-1.15302,true,1),
        		new ScoredPSM(-1.15922,true,1),
        		new ScoredPSM(-1.16116,true,1),
        		new ScoredPSM(-1.16129,true,1),
        		new ScoredPSM(-1.16349,true,1),
        		new ScoredPSM(-1.16513,true,1),
        		new ScoredPSM(-1.16519,true,1),
        		new ScoredPSM(-1.16648,true,1),
        		new ScoredPSM(-1.16833,true,1),
        		new ScoredPSM(-1.16866,true,1),
        		new ScoredPSM(-1.1705,true,1),
        		new ScoredPSM(-1.17133,true,1),
        		new ScoredPSM(-1.17207,true,1),
        		new ScoredPSM(-1.174,true,1),
        		new ScoredPSM(-1.17653,true,1),
        		new ScoredPSM(-1.17702,true,1),
        		new ScoredPSM(-1.17775,true,1),
        		new ScoredPSM(-1.1781,true,1),
        		new ScoredPSM(-1.17966,true,1),
        		new ScoredPSM(-1.17968,true,1),
        		new ScoredPSM(-1.17996,true,1),
        		new ScoredPSM(-1.18668,true,1),
        		new ScoredPSM(-1.19201,true,1),
        		new ScoredPSM(-1.1931,true,1),
        		new ScoredPSM(-1.19648,true,1),
        		new ScoredPSM(-1.19838,true,1),
        		new ScoredPSM(-1.19901,true,1),
        		new ScoredPSM(-1.20208,true,1),
        		new ScoredPSM(-1.20317,true,1),
        		new ScoredPSM(-1.20368,true,1),
        		new ScoredPSM(-1.20464,true,1),
        		new ScoredPSM(-1.20861,true,1),
        		new ScoredPSM(-1.21249,true,1),
        		new ScoredPSM(-1.21522,true,1),
        		new ScoredPSM(-1.21867,true,1),
        		new ScoredPSM(-1.22247,true,1),
        		new ScoredPSM(-1.22542,true,1),
        		new ScoredPSM(-1.22574,true,1),
        		new ScoredPSM(-1.22939,true,1),
        		new ScoredPSM(-1.23119,true,1),
        		new ScoredPSM(-1.23577,true,1),
        		new ScoredPSM(-1.23862,true,1),
        		new ScoredPSM(-1.23903,true,1),
        		new ScoredPSM(-1.24145,true,1),
        		new ScoredPSM(-1.25204,true,1),
        		new ScoredPSM(-1.25902,true,1),
        		new ScoredPSM(-1.26033,true,1),
        		new ScoredPSM(-1.26208,true,1),
        		new ScoredPSM(-1.26258,true,1),
        		new ScoredPSM(-1.26407,true,1),
        		new ScoredPSM(-1.26854,true,1),
        		new ScoredPSM(-1.26952,true,1),
        		new ScoredPSM(-1.27119,true,1),
        		new ScoredPSM(-1.27414,true,1),
        		new ScoredPSM(-1.27575,true,1),
        		new ScoredPSM(-1.27626,true,1),
        		new ScoredPSM(-1.27891,true,1),
        		new ScoredPSM(-1.28058,true,1),
        		new ScoredPSM(-1.28264,true,1),
        		new ScoredPSM(-1.28345,true,1),
        		new ScoredPSM(-1.28349,true,1),
        		new ScoredPSM(-1.28371,true,1),
        		new ScoredPSM(-1.28471,true,1),
        		new ScoredPSM(-1.29588,true,1),
        		new ScoredPSM(-1.30058,true,1),
        		new ScoredPSM(-1.30087,true,1),
        		new ScoredPSM(-1.30234,true,1),
        		new ScoredPSM(-1.30431,true,1),
        		new ScoredPSM(-1.30589,true,1),
        		new ScoredPSM(-1.30624,true,1),
        		new ScoredPSM(-1.30855,true,1),
        		new ScoredPSM(-1.31167,true,1),
        		new ScoredPSM(-1.31663,true,1),
        		new ScoredPSM(-1.31696,true,1),
        		new ScoredPSM(-1.31747,true,1),
        		new ScoredPSM(-1.31835,true,1),
        		new ScoredPSM(-1.32254,true,1),
        		new ScoredPSM(-1.32514,true,1),
        		new ScoredPSM(-1.32637,true,1),
        		new ScoredPSM(-1.33535,true,1),
        		new ScoredPSM(-1.3355,true,1),
        		new ScoredPSM(-1.34076,true,1),
        		new ScoredPSM(-1.34899,true,1),
        		new ScoredPSM(-1.35022,true,1),
        		new ScoredPSM(-1.35526,true,1),
        		new ScoredPSM(-1.36355,true,1),
        		new ScoredPSM(-1.36974,true,1),
        		new ScoredPSM(-1.36985,true,1),
        		new ScoredPSM(-1.37059,true,1),
        		new ScoredPSM(-1.37645,true,1),
        		new ScoredPSM(-1.38141,true,1),
        		new ScoredPSM(-1.38611,true,1),
        		new ScoredPSM(-1.38893,true,1),
        		new ScoredPSM(-1.39466,true,1),
        		new ScoredPSM(-1.39669,true,1),
        		new ScoredPSM(-1.39738,true,1),
        		new ScoredPSM(-1.40314,true,1),
        		new ScoredPSM(-1.40808,true,1),
        		new ScoredPSM(-1.41406,true,1),
        		new ScoredPSM(-1.41473,true,1),
        		new ScoredPSM(-1.41796,true,1),
        		new ScoredPSM(-1.41933,true,1),
        		new ScoredPSM(-1.42133,true,1),
        		new ScoredPSM(-1.42228,true,1),
        		new ScoredPSM(-1.42696,true,1),
        		new ScoredPSM(-1.42794,true,1),
        		new ScoredPSM(-1.42906,true,1),
        		new ScoredPSM(-1.43438,true,1),
        		new ScoredPSM(-1.43507,true,1),
        		new ScoredPSM(-1.43542,true,1),
        		new ScoredPSM(-1.44042,true,1),
        		new ScoredPSM(-1.44266,true,1),
        		new ScoredPSM(-1.44922,true,1),
        		new ScoredPSM(-1.45401,true,1),
        		new ScoredPSM(-1.45513,true,1),
        		new ScoredPSM(-1.4571,true,1),
        		new ScoredPSM(-1.45825,true,1),
        		new ScoredPSM(-1.4586,true,1),
        		new ScoredPSM(-1.46696,true,1),
        		new ScoredPSM(-1.46709,true,1),
        		new ScoredPSM(-1.46744,true,1),
        		new ScoredPSM(-1.4787,true,1),
        		new ScoredPSM(-1.48283,true,1),
        		new ScoredPSM(-1.48482,true,1),
        		new ScoredPSM(-1.48494,true,1),
        		new ScoredPSM(-1.48874,true,1),
        		new ScoredPSM(-1.48904,true,1),
        		new ScoredPSM(-1.48974,true,1),
        		new ScoredPSM(-1.4898,true,1),
        		new ScoredPSM(-1.4929,true,1),
        		new ScoredPSM(-1.49496,true,1),
        		new ScoredPSM(-1.49857,true,1),
        		new ScoredPSM(-1.49889,true,1),
        		new ScoredPSM(-1.50176,true,1),
        		new ScoredPSM(-1.50653,true,1),
        		new ScoredPSM(-1.50676,true,1),
        		new ScoredPSM(-1.50747,true,1),
        		new ScoredPSM(-1.51309,true,1),
        		new ScoredPSM(-1.51443,true,1),
        		new ScoredPSM(-1.51882,true,1),
        		new ScoredPSM(-1.52363,true,1),
        		new ScoredPSM(-1.52922,true,1),
        		new ScoredPSM(-1.52957,true,1),
        		new ScoredPSM(-1.53388,true,1),
        		new ScoredPSM(-1.53455,true,1),
        		new ScoredPSM(-1.53564,true,1),
        		new ScoredPSM(-1.54653,true,1),
        		new ScoredPSM(-1.54837,true,1),
        		new ScoredPSM(-1.56153,true,1),
        		new ScoredPSM(-1.56219,true,1),
        		new ScoredPSM(-1.56438,true,1),
        		new ScoredPSM(-1.57289,true,1),
        		new ScoredPSM(-1.58373,true,1),
        		new ScoredPSM(-1.58728,true,1),
        		new ScoredPSM(-1.59375,true,1),
        		new ScoredPSM(-1.60495,true,1),
        		new ScoredPSM(-1.60583,true,1),
        		new ScoredPSM(-1.6186,true,1),
        		new ScoredPSM(-1.62339,true,1),
        		new ScoredPSM(-1.63349,true,1),
        		new ScoredPSM(-1.64018,true,1),
        		new ScoredPSM(-1.64891,true,1),
        		new ScoredPSM(-1.65164,true,1),
        		new ScoredPSM(-1.65362,true,1),
        		new ScoredPSM(-1.65633,true,1),
        		new ScoredPSM(-1.66937,true,1),
        		new ScoredPSM(-1.66947,true,1),
        		new ScoredPSM(-1.68084,true,1),
        		new ScoredPSM(-1.68206,true,1),
        		new ScoredPSM(-1.6822,true,1),
        		new ScoredPSM(-1.69064,true,1),
        		new ScoredPSM(-1.69541,true,1),
        		new ScoredPSM(-1.69851,true,1),
        		new ScoredPSM(-1.70105,true,1),
        		new ScoredPSM(-1.70583,true,1),
        		new ScoredPSM(-1.71253,true,1),
        		new ScoredPSM(-1.73277,true,1),
        		new ScoredPSM(-1.76169,true,1),
        		new ScoredPSM(-1.76554,true,1),
        		new ScoredPSM(-1.76581,true,1),
        		new ScoredPSM(-1.76702,true,1),
        		new ScoredPSM(-1.76972,true,1),
        		new ScoredPSM(-1.77645,true,1),
        		new ScoredPSM(-1.7772,true,1),
        		new ScoredPSM(-1.7792,true,1),
        		new ScoredPSM(-1.78219,true,1),
        		new ScoredPSM(-1.78394,true,1),
        		new ScoredPSM(-1.82344,true,1),
        		new ScoredPSM(-1.82655,true,1),
        		new ScoredPSM(-1.83204,true,1),
        		new ScoredPSM(-1.8327,true,1),
        		new ScoredPSM(-1.84458,true,1),
        		new ScoredPSM(-1.84627,true,1),
        		new ScoredPSM(-1.87223,true,1),
        		new ScoredPSM(-1.87925,true,1),
        		new ScoredPSM(-1.89082,true,1),
        		new ScoredPSM(-1.9071,true,1),
        		new ScoredPSM(-1.91069,true,1),
        		new ScoredPSM(-1.92433,true,1),
        		new ScoredPSM(-1.93662,true,1),
        		new ScoredPSM(-1.94326,true,1),
        		new ScoredPSM(-1.96667,true,1),
        		new ScoredPSM(-2.03616,true,1),
        		new ScoredPSM(-2.04246,true,1),
        		new ScoredPSM(-2.04732,true,1),
        		new ScoredPSM(-2.06739,true,1),
        		new ScoredPSM(-2.12761,true,1),
        		new ScoredPSM(-2.1343,true,1),
        		new ScoredPSM(-2.17002,true,1),
        		new ScoredPSM(-2.18661,true,1),
        		new ScoredPSM(-2.2337,true,1),
        		new ScoredPSM(-2.29051,true,1),
        		new ScoredPSM(-2.31747,true,1),
        		new ScoredPSM(-2.32941,true,1),
        		new ScoredPSM(-2.33445,true,1),
        		new ScoredPSM(-2.33729,true,1),
        		new ScoredPSM(-2.43982,true,1),
        		new ScoredPSM(-2.44553,true,1),
        		new ScoredPSM(-2.45925,true,1),
        		new ScoredPSM(-2.50692,true,1),
        		new ScoredPSM(-2.53972,true,1),
        		new ScoredPSM(-2.54228,true,1),
        		new ScoredPSM(-2.54578,true,1),
        		new ScoredPSM(-2.5484,true,1),
        		new ScoredPSM(-2.55152,true,1),
        		new ScoredPSM(-2.57304,true,1),
        		new ScoredPSM(-2.68186,true,1),
        		new ScoredPSM(-2.69152,true,1),
        		new ScoredPSM(-2.77612,true,1),
        		new ScoredPSM(-2.77736,true,1)
        );
        Collections.sort(psms);
        
        Function function = d2pep(psms);
        
        int numberOfTargets=0;
        int numberOfDecoys=0;
        double estimatedTruesOriginally=0.0;
        double estimatedTrues=0.0;
        for (ScoredPSM scoredPSM : psms) {
        	if (!scoredPSM.isDecoy()) {
        		estimatedTruesOriginally+=1.0-scoredPSM.getEstimatedPEP();
        		estimatedTrues+=1.0-function.getYValue((float)scoredPSM.getPrimaryScore());
        		numberOfTargets++;
        	} else {
        		numberOfDecoys++;
        	}
		}
        System.out.println("T:"+numberOfTargets+", D:"+numberOfDecoys+", estimated Trues:"+estimatedTrues+", estimated Trues Originally:"+estimatedTruesOriginally);
        
        //System.out.println("PEP (targets+decoys, input order): " + peps);
    }
}
