package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

public class LaplaceDistributionTest {
	public static void main(String[] args) {
		Distribution d=new LaplaceDistribution(-0.5, 1, 1);
		for (double x = -5; x < 5; x+=0.1) {
			System.out.println(x+"\t"+d.getProbability(x)+"\t"+d.getCDF(x));
		}
	}
}
