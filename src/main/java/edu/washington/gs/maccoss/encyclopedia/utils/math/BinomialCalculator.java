package edu.washington.gs.maccoss.encyclopedia.utils.math;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.special.Beta;

public class BinomialCalculator extends BinomialDistribution {
	private static final long serialVersionUID=1L;

	public BinomialCalculator(int trials, double p) {
		super(new Well19937c(), trials, p);
	}

    /** {@inheritDoc} */
    public double cumulativeProbabilityGreaterThan(int x) {
        double ret;
        if (x < 0) {
            ret = 0.0;
        } else if (x >= getNumberOfTrials()) {
            ret = 1.0;
        } else {
            ret = Beta.regularizedBeta(getProbabilityOfSuccess(), x + 1.0, getNumberOfTrials() - x);
        }
        return ret;
    }
}
