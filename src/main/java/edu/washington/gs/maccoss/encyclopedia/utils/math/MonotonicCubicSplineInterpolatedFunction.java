package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class MonotonicCubicSplineInterpolatedFunction implements Function {
	private final ArrayList<XYPoint> knots;
	private final double[] x;
	private final double[] yIso;   // isotonic fit of y (increasing or decreasing)
	private final double[] h;      // segment widths
	private final double[] delta;  // segment secant slopes
	private final double[] m;      // node derivatives (Fritsch–Carlson)
	private final boolean increasing;

	public MonotonicCubicSplineInterpolatedFunction(ArrayList<XYPoint> knots, boolean increasing) {
		this.knots = knots;
		this.increasing = increasing;

		Pair<double[], double[]> xys = XYTrace.toArrays(knots);
		this.x = xys.x;

		// Enforce chosen monotonicity using a sign trick:
		// if decreasing, negate y, run PAVA non-decreasing, then flip back.
		double sign = increasing ? 1.0 : -1.0;
		double[] yRaw = Arrays.copyOf(xys.y, xys.y.length);
		for (int i = 0; i < yRaw.length; i++) yRaw[i] *= sign;

		pavaNonDecreasingInPlace(yRaw);  // now yRaw is non-decreasing
		this.yIso = new double[yRaw.length];
		for (int i = 0; i < yRaw.length; i++) this.yIso[i] = sign * yRaw[i];

		// Precompute segments and monotone Hermite derivatives
		final int n = x.length;
		this.h = new double[Math.max(0, n - 1)];
		this.delta = new double[Math.max(0, n - 1)];
		this.m = new double[Math.max(0, n)];

		if (n >= 2) {
			for (int i = 0; i < n - 1; i++) {
				h[i] = x[i + 1] - x[i];
				if (h[i] <= 0) {
					// Degenerate or unsorted x; guard to avoid NaNs
					h[i] = 1e-12;
				}
				delta[i] = (yIso[i + 1] - yIso[i]) / h[i];
			}

			// Endpoint slopes
			m[0] = delta[0];
			m[n - 1] = delta[n - 2];

			// Interior slopes with Fritsch–Carlson formula
			for (int i = 1; i < n - 1; i++) {
				final double d1 = delta[i - 1];
				final double d2 = delta[i];
				if (d1 * d2 <= 0.0) {
					m[i] = 0.0;
				} else {
					double w1 = 2.0 * h[i] + h[i - 1];
					double w2 = h[i] + 2.0 * h[i - 1];
					m[i] = (w1 + w2) / (w1 / d1 + w2 / d2);
				}
			}

			// Extra safety: if a segment is flat, zero the adjacent slopes
			for (int i = 0; i < n - 1; i++) {
				if (Math.abs(delta[i]) == 0.0) {
					m[i] = 0.0;
					m[i + 1] = 0.0;
				}
			}
		}
	}

	@Override
	public ArrayList<XYPoint> getKnots() {
		return knots;
	}

	@Override
	public boolean isXInsideBoundaries(float xi) {
		int upperBin = calculateBinNumberIncreasing(xi, x);
		if (upperBin == 0) return false;
		if (upperBin == x.length) return false;
		return true;
	}

	@Override
	public float getYValue(float xi) {
		int upperBin = calculateBinNumberIncreasing(xi, x);

		// boundary conditions
		if (upperBin <= 0) return (float) yIso[0];
		if (upperBin >= x.length) return (float) yIso[yIso.length - 1];

		int i = upperBin - 1;
		double xl = x[i];
		double xr = x[i + 1];
		double t = (xi - xl) / (xr - xl);
		if (t <= 0) return (float) yIso[i];
		if (t >= 1) return (float) yIso[i + 1];

		return (float) cubicHermite(yIso[i], yIso[i + 1], m[i], m[i + 1], xr - xl, t);
	}

	@Override
	public boolean isYInsideBoundaries(float yi) {
		// Works for both directions
		double yMin = Math.min(yIso[0], yIso[yIso.length - 1]);
		double yMax = Math.max(yIso[0], yIso[yIso.length - 1]);
		return yi > yMin && yi < yMax;
	}

	@Override
	public float getXValue(float yi) {
		// boundary conditions
		if (increasing) {
			if (yi <= yIso[0]) return (float) x[0];
			if (yi >= yIso[yIso.length - 1]) return (float) x[x.length - 1];
		} else {
			if (yi >= yIso[0]) return (float) x[0];
			if (yi <= yIso[yIso.length - 1]) return (float) x[x.length - 1];
		}

		// Find the knot interval in y-space respecting monotonic direction
		int upperBin = calculateBinNumberMonotone(yi, yIso, increasing);
		int i = Math.max(0, Math.min(upperBin - 1, x.length - 2));
		double xl = x[i];
		double xr = x[i + 1];

		// If the segment is exactly flat, return midpoint
		if (yIso[i] == yIso[i + 1]) {
			return (float) ((xl + xr) * 0.5);
		}

		// Invert F(x)=yi on [xl, xr] by monotone bisection using the Hermite segment
		double left = xl, right = xr;
		for (int iter = 0; iter < 40; iter++) { // ~1e-12 relative with well-scaled x
			double mid = 0.5 * (left + right);
			double t = (mid - xl) / (xr - xl);
			double fmid = cubicHermite(yIso[i], yIso[i + 1], m[i], m[i + 1], xr - xl, t);

			// Direction-aware comparison
			if ((increasing && fmid < yi) || (!increasing && fmid > yi)) {
				left = mid;
			} else {
				right = mid;
			}
		}
		return (float) (0.5 * (left + right));
	}

	/* ===== Helpers ===== */

	/**
	 * Cubic Hermite (shape-preserving) evaluation on one segment.
	 * y(t) = H00*y0 + H10*m0*h + H01*y1 + H11*m1*h, with t in [0,1]
	 */
	private static double cubicHermite(double y0, double y1, double m0, double m1, double h, double t) {
		double t2 = t * t;
		double t3 = t2 * t;
		double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
		double h10 = t3 - 2.0 * t2 + t;
		double h01 = -2.0 * t3 + 3.0 * t2;
		double h11 = t3 - t2;
		return h00 * y0 + h10 * m0 * h + h01 * y1 + h11 * m1 * h;
	}

	/**
	 * Pool-Adjacent-Violators Algorithm: enforce y non-decreasing in-place.
	 * (All weights = 1 here.)
	 */
	private static void pavaNonDecreasingInPlace(double[] y) {
		final int n = y.length;
		if (n <= 1) return;

		double[] v = Arrays.copyOf(y, n);
		int[]    c = new int[n]; // counts (weights as integers)
		int m = -1;              // last active block index

		for (int i = 0; i < n; i++) {
			// start a new block with value y[i], count 1
			m++;
			v[m] = y[i];
			c[m] = 1;

			// merge backward while violating non-decreasing order
			while (m > 0 && v[m - 1] > v[m]) {
				int totC = c[m - 1] + c[m];
				v[m - 1] = (c[m - 1] * v[m - 1] + c[m] * v[m]) / totC;
				c[m - 1] = totC;
				m--;
			}
		}

		// expand block averages back to y
		int idx = 0;
		for (int b = 0; b <= m; b++) {
			for (int k = 0; k < c[b]; k++) {
				y[idx++] = v[b];
			}
		}
	}

	/**
	 * Binary search for insertion point assuming xs is **increasing**.
	 * Returns the same-style "upper bin" index used in your code.
	 */
	public static int calculateBinNumberIncreasing(double x, double[] xs) {
		int bin = Arrays.binarySearch(xs, x);
		if (bin < 0) bin = (-(bin + 1));
		return bin;
	}

	/**
	 * Binary search for insertion point in a monotone (increasing OR decreasing) array.
	 * Returns the "upper bin" index consistent with your usage.
	 */
	public static int calculateBinNumberMonotone(double value, double[] arr, boolean increasing) {
		if (increasing) {
			return calculateBinNumberIncreasing(value, arr);
		} else {
			// Manual binary search for decreasing array
			int lo = 0, hi = arr.length; // insertion point in [0..n]
			while (lo < hi) {
				int mid = (lo + hi) >>> 1;
				// In decreasing order: arr[mid] >= value means move right boundary left
				if (arr[mid] >= value) {
					lo = mid + 1;
				} else {
					hi = mid;
				}
			}
			// For decreasing arrays, lo is the count of elements >= value, i.e., the "upper bin"
			return lo;
		}
	}
}
