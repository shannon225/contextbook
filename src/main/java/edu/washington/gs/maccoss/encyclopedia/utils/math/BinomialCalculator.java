package edu.washington.gs.maccoss.encyclopedia.utils.math;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.special.Beta;

public class BinomialCalculator extends BinomialDistribution {
	private static final long serialVersionUID=1L;
	
	public static void main(String[] args) {
		System.out.println("Ascore 1 tolerance, peakdepth 10");
		BinomialCalculator calc=new BinomialCalculator(8, 0.1);
		double correct = calc.cumulativeProbabilityGreaterThan(4);
		double incorrect = calc.cumulativeProbabilityGreaterThan(0);
		System.out.println(correct+", "+incorrect+", "+-10*(Log.log10(correct)-Log.log10(incorrect)));
		System.out.println();

		System.out.println("Ascore 1 tolerance, peakdepth 1");
		calc=new BinomialCalculator(8, 0.01);
		correct = calc.cumulativeProbabilityGreaterThan(1);
		incorrect = 1.0f;
		System.out.println(correct+", "+incorrect+", "+-10*(Log.log10(correct)-Log.log10(incorrect)));
		System.out.println();

		System.out.println("Ascore 0.1 tolerance, peakdepth 10");
		calc=new BinomialCalculator(8, 0.01);
		correct = calc.cumulativeProbabilityGreaterThan(4);
		incorrect = calc.cumulativeProbabilityGreaterThan(0);
		System.out.println(correct+", "+incorrect+", "+-10*(Log.log10(correct)-Log.log10(incorrect)));
		System.out.println();
		
		System.out.println("Peptide Score");
		calc=new BinomialCalculator(18, 0.1);
		correct = calc.cumulativeProbabilityGreaterThan(8);
		incorrect = calc.cumulativeProbabilityGreaterThan(5);
		System.out.println(correct+", "+incorrect+", "+-10*(Log.log10(correct)-Log.log10(incorrect)));
		System.out.println();
	}

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
