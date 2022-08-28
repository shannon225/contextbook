package edu.washington.gs.maccoss.encyclopedia.utils.math;

/**
 * based on gradient descent of the algebraic sigmoid described here:
 * http://rstudio-pubs-static.s3.amazonaws.com/252141_6c6b1b2857a04a80a6bbfbc0481e1469.html
 * Formulation is (a+b*x)*c/sqrt(1+(a+b*x)^2) + d
 * 
 * This is better than a logistic function for fitting when the min and max aren't known (e.g., fitting pi0)
 * @author searleb
 *
 */
public class Sigmoid {
	protected float a;
	protected float b;
	protected float c;
	protected float d;
	private float alpha;
	
	public Sigmoid() {
		this(0.0f, 0.0f, 0.0f, 0.0f);
	}

	public Sigmoid(float alpha) {
		this(0.0f, 0.0f, 0.0f, 0.0f, alpha);
	}
	
	public Sigmoid(float a, float b, float c, float d) {
		this(a, b, c, d, 0.0000001f);
	}

	public Sigmoid(float a, float b, float c, float d, float alpha) {
		this.a = a;
		this.b = b;
		this.c = c;
		this.d = d;
		this.alpha = alpha;
	}
	
	public void train(float[] xs, float[] ys, int maxIterations) {
		for (int i = 0; i < maxIterations; i++) {
			float oldA=a;
			float oldB=b;
			float oldC=c;
			float oldD=d;
			
			float gradientA = alpha*gradientA(xs, ys);
			adjustA(gradientA);
			float gradientB = alpha*gradientB(xs, ys);
			adjustB(gradientB);
			float gradientC = alpha*gradientC(xs, ys);
			adjustC(gradientC);
			float gradientD = alpha*gradientD(xs, ys);
			adjustD(gradientD);
			
			// change beyond floating point error
			if (a==oldA&&b==oldB&&c==oldC&&d==oldD) break;
		}
	}

	protected float adjustD(float gradientD) {
		return d=d-gradientD;
	}

	protected float adjustC(float gradientC) {
		return c=c-gradientC;
	}

	protected float adjustB(float gradientB) {
		return b=b-gradientB;
	}

	protected float adjustA(float gradientA) {
		return a=a-gradientA;
	}
	
	/**
	 * sum(   -2*c*( (y - d)*sqrt((a+b*x)^2 +1 )-c*(a+b*x)) /   ((a+b*x)^2 + 1)^2  )
	 * @param xs
	 * @param ys
	 * @return
	 */
	private float gradientA(float[] xs, float[] ys) {
		float[] xAdjs=General.add(General.multiply(xs, b), a);
		float[] xAdjsSqaredPlus1=General.add(General.multiply(xAdjs, xAdjs), 1.0f);
		float[] numerator=General.subtract(General.multiply(General.subtract(ys, d), General.protectedSqrt(xAdjsSqaredPlus1)), General.multiply(xAdjs, c));
		return General.sum(General.divide(General.multiply(numerator, -2.0f*c), General.multiply(xAdjsSqaredPlus1, xAdjsSqaredPlus1)));
	}
	
	/**
	 * sum(   -2*c*x*( (y - d)*sqrt((a+b*x)^2 +1 )-c*(a+b*x)) /   ((a+b*x)^2 + 1)^2  )
	 * @param xs
	 * @param ys
	 * @return
	 */
	private float gradientB(float[] xs, float[] ys) {
		float[] xAdjs=General.add(General.multiply(xs, b), a);
		float[] xAdjsSqaredPlus1=General.add(General.multiply(xAdjs, xAdjs), 1.0f);
		float[] numerator=General.subtract(General.multiply(General.subtract(ys, d), General.protectedSqrt(xAdjsSqaredPlus1)), General.multiply(xAdjs, c));
		return General.sum(General.divide(General.multiply(numerator, General.multiply(xs, -2.0f*c)), General.multiply(xAdjsSqaredPlus1, xAdjsSqaredPlus1)));
	}

	/**
	 * sum(  (2*((b*x+a)*c- sqrt((b*x+a)^2+1)*(y-d)))/((b*x+a)+1)   )
	 * @param xs
	 * @param ys
	 * @return
	 */
	private float gradientC(float[] xs, float[] ys) {
		float[] xAdjs=General.add(General.multiply(xs, b), a);
		float[] xAdjsSqaredPlus1=General.add(General.multiply(xAdjs, xAdjs), 1.0f);
		
		float[] numerator=General.subtract(General.multiply(xAdjs, c), General.multiply(General.protectedSqrt(xAdjsSqaredPlus1), General.subtract(ys,  d)));
		return General.sum(General.divide(General.multiply(numerator, 2.0f), General.add(xAdjs, 1.0f)));
	}

	/**
	 * sum( -2*(-d+y -(c*(b*x+a))/sqrt((b*x+a)^2 + 1)    ) )
	 * @param xs
	 * @param ys
	 * @return
	 */
	private float gradientD(float[] xs, float[] ys) {
		float[] xAdjs=General.add(General.multiply(xs, b), a);
		float[] xAdjsSqaredPlus1=General.add(General.multiply(xAdjs, xAdjs), 1.0f);
		
		float[] subTerm=General.divide(General.multiply(xAdjs, c), General.protectedSqrt(xAdjsSqaredPlus1));
		return General.sum(General.multiply(General.subtract(General.subtract(ys, d), subTerm), -2.0f));
	}
	
	public float getValue(float x) {
		float xAdj=a+b*x;
		return xAdj*c/(float)Math.sqrt(1+xAdj*xAdj) + d;
	}

	public float[] getValues(float[] xs) {
		float[] ys=new float[xs.length];
		for (int i = 0; i < ys.length; i++) {
			ys[i]=getValue(xs[i]);
		}
		return ys;
	}
	
	public float getA() {
		return a;
	}
	
	public float getB() {
		return b;
	}
	
	public float getC() {
		return c;
	}
	
	public float getD() {
		return d;
	}
}
