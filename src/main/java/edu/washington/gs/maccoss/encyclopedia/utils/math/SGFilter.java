package edu.washington.gs.maccoss.encyclopedia.utils.math;

/*
 * Heavily adapted from:
 * https://github.com/swallez/savitzky-golay-filter/tree/master/src/mr/go/sgfilter
 * 
 * Copyright [2009] [Marcin Rzeźnicki]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
 */

import java.util.Arrays;

//import org.apache.commons.math.linear.RealMatrixImpl
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealVector;

/**
 * Savitzky-Golay filter implementation. For more information see
 * http://www.nrbook.com/a/bookcpdf/c14-8.pdf. This implementation, however,
 * does not use FFT Note: modified version by DBO for compatibility with
 * apache.commons.math3
 * 
 * @author Marcin Rzeźnicki
 * 
 */
public class SGFilter {

	/**
	 * Computes Savitzky-Golay coefficients for given parameters
	 * 
	 * @param nl
	 *            numer of past data points filter will use
	 * @param nr
	 *            number of future data points filter will use
	 * @param degree
	 *            order of smoothin polynomial
	 * @return Savitzky-Golay coefficients
	 * @throws IllegalArgumentException
	 *             if {@code nl < 0} or {@code nr < 0} or {@code nl + nr <
	 *             degree}
	 */
	private static float[] computeSGCoefficients(int nl, int nr, int degree) {
		// RealMatrixImpl matrix = new RealMatrixImpl(degree + 1, degree + 1);
		Array2DRowRealMatrix matrix=new Array2DRowRealMatrix(degree+1, degree+1);
		double[][] a=matrix.getDataRef();
		float sum;
		for (int i=0; i<=degree; i++) {
			for (int j=0; j<=degree; j++) {
				sum=(i==0&&j==0)?1:0;
				for (int k=1; k<=nr; k++)
					sum+=Math.pow(k, i+j);
				for (int k=1; k<=nl; k++)
					sum+=Math.pow(-k, i+j);
				a[i][j]=sum;
			}
		}

		double[] b=new double[degree+1];
		b[0]=1;

		final RealVector constantVector=new ArrayRealVector(b, false);
		DecompositionSolver solver=new LUDecomposition(matrix).getSolver();

		b=solver.solve(constantVector).toArray();

		// b = matrix.solve(b);
		float[] coeffs=new float[nl+nr+1];
		for (int n=-nl; n<=nr; n++) {
			sum=(float)b[0];
			for (int m=1; m<=degree; m++)
				sum+=b[m]*Math.pow(n, m);
			coeffs[n+nl]=sum;
		}
		return coeffs;
	}

	private final int nl;
	private final int nr;
	private final float[][] coefficients;

	/**
	 * Constructs Savitzky-Golay filter which uses specified numebr of
	 * surrounding data points
	 * 
	 * @param nl
	 *            numer of past data points filter will use
	 * @param nr
	 *            numer of future data points filter will use
	 * @param degree
	 *            order of smoothin polynomial
	 * @throws IllegalArgumentException
	 *             of {@code nl < 0} or {@code nr < 0}
	 */
	public SGFilter(int nl, int nr, int degree) {
		if (nl<0||nr<0||nl+nr<degree) throw new IllegalArgumentException("Bad arguments");
		this.nl=nl;
		this.nr=nr;
		this.coefficients=new float[][] {computeSGCoefficients(nl, nr, degree)};
	}

	/**
	 * Smooths data by using Savitzky-Golay filter. This method will use 0 for
	 * any element beyond {@code data} which will be needed for computation
	 * 
	 * @param data
	 *            data for filter
	 * @param coeffs
	 *            filter coefficients
	 * @return filtered data
	 * @throws NullPointerException
	 *             when any array passed as parameter is null
	 */
	public float[] smooth(float[] data) {
		return smooth(data, 0, data.length);
	}

	/**
	 * Smooths data by using Savitzky-Golay filter. Smoothing uses {@code
	 * leftPad} and/or {@code rightPad} if you want to augment data on
	 * boundaries to achieve smoother results for your purpose. If you do not
	 * need this feature you may pass empty arrays (filter will use 0s in this
	 * place, so you may want to use appropriate preprocessor)
	 * 
	 * @param data
	 *            data for filter
	 * @param leftPad
	 *            left padding
	 * @param rightPad
	 *            right padding
	 * @param coeffs
	 *            filter coefficients
	 * @return filtered data
	 * @throws NullPointerException
	 *             when any array passed as parameter is null
	 */
	public float[] smooth(float[] data, float[] leftPad, float[] rightPad) {
		return smooth(data, leftPad, rightPad, 0);
	}

	/**
	 * Smooths data by using Savitzky-Golay filter. Smoothing uses {@code
	 * leftPad} and/or {@code rightPad} if you want to augment data on
	 * boundaries to achieve smoother results for your purpose. If you do not
	 * need this feature you may pass empty arrays (filter will use 0s in this
	 * place, so you may want to use appropriate preprocessor). If you want to
	 * use different (probably non-symmetrical) filter near both ends of
	 * (padded) data, you will be using {@code bias} and {@code coeffs}. {@code
	 * bias} essentially means
	 * "how many points of pad should be left out when smoothing". Filters
	 * taking this condition into consideration are passed in {@code coeffs}.
	 * <tt>coeffs[0]</tt> is used for unbiased data (that is, for
	 * <tt>data[bias]..data[data.length-bias-1]</tt>). Its length has to be
	 * <tt>nr + nl + 1</tt>. Filters from range
	 * <tt>coeffs[coeffs.length - 1]</tt> to
	 * <tt>coeffs[coeffs.length - bias]</tt> are used for smoothing first
	 * {@code bias} points (that is, from <tt>data[0]</tt> to
	 * <tt>data[bias]</tt>) correspondingly. Filters from range
	 * <tt>coeffs[1]</tt> to <tt>coeffs[bias]</tt> are used for smoothing last
	 * {@code bias} points (that is, for
	 * <tt>data[data.length-bias]..data[data.length-1]</tt>). For example, if
	 * you use 5 past points and 5 future points for smoothing, but have only 3
	 * meaningful padding points - you would use {@code bias} equal to 2 and
	 * would pass in {@code coeffs} param filters taking 5-5 points (for regular
	 * smoothing), 5-4, 5-3 (for rightmost range of data) and 3-5, 4-5 (for
	 * leftmost range). If you do not wish to use pads completely for
	 * symmetrical filter then you should pass <tt>bias = nl = nr</tt>
	 * 
	 * @param data
	 *            data for filter
	 * @param leftPad
	 *            left padding
	 * @param rightPad
	 *            right padding
	 * @param bias
	 *            how many points of pad should be left out when smoothing
	 * @param coeffs
	 *            array of filter coefficients
	 * @return filtered data
	 * @throws IllegalArgumentException
	 *             when <tt>bias < 0</tt> or <tt>bias > min(nr, nl)</tt>
	 * @throws IndexOutOfBoundsException
	 *             when {@code coeffs} has less than <tt>2*bias + 1</tt>
	 *             elements
	 * @throws NullPointerException
	 *             when any array passed as parameter is null
	 */
	public float[] smooth(float[] data, float[] leftPad, float[] rightPad, int bias) {
		if (bias<0||bias>nr||bias>nl) throw new IllegalArgumentException("bias < 0 or bias > nr or bias > nl");
		int dataLength=data.length;
		if (dataLength==0) return data;
		int n=dataLength+nl+nr;
		float[] dataCopy=new float[n];
		// copy left pad reversed
		int leftPadOffset=nl-leftPad.length;
		if (leftPadOffset>=0) for (int i=0; i<leftPad.length; i++) {
			dataCopy[leftPadOffset+i]=leftPad[i];
		}
		else for (int i=0; i<nl; i++) {
			dataCopy[i]=leftPad[i-leftPadOffset];
		}
		// copy actual data
		for (int i=0; i<dataLength; i++) {
			dataCopy[i+nl]=data[i];
		}
		// copy right pad
		int rightPadOffset=nr-rightPad.length;
		if (rightPadOffset>=0) for (int i=0; i<rightPad.length; i++) {
			dataCopy[i+dataLength+nl]=rightPad[i];
		}
		else for (int i=0; i<nr; i++) {
			dataCopy[i+dataLength+nl]=rightPad[i];
		}

		// convolution (with savitzky-golay coefficients)
		float[] sdata=new float[dataLength];
		float[] sg;
		for (int b=bias; b>0; b--) {
			sg=coefficients[coefficients.length-b];
			int x=(nl+bias)-b;
			float sum=0;
			for (int i=-nl+b; i<=nr; i++) {
				sum+=dataCopy[x+i]*sg[nl-b+i];
			}
			sdata[x-nl]=sum;
		}
		sg=coefficients[0];
		for (int x=nl+bias; x<n-nr-bias; x++) {
			float sum=0;
			for (int i=-nl; i<=nr; i++) {
				sum+=dataCopy[x+i]*sg[nl+i];
			}
			sdata[x-nl]=sum;
		}
		for (int b=1; b<=bias; b++) {
			sg=coefficients[b];
			int x=(n-nr-bias)+(b-1);
			float sum=0;
			for (int i=-nl; i<=nr-b; i++) {
				sum+=dataCopy[x+i]*sg[nl+i];
			}
			sdata[x-nl]=sum;
		}
		return sdata;
	}

	/**
	 * Runs filter on data from {@code from} (including) to {@code to}
	 * (excluding). Data beyond range spanned by {@code from} and {@code to}
	 * will be used for padding
	 * 
	 * @param data
	 *            data for filter
	 * @param from
	 *            inedx of the first element of data
	 * @param to
	 *            index of the first element omitted
	 * @param coeffs
	 *            filter coefficients
	 * @return filtered data
	 * @throws ArrayIndexOutOfBoundsException
	 *             if <tt>to > data.length</tt>
	 * @throws IllegalArgumentException
	 *             if <tt>from < 0</tt> or <tt>to > data.length</tt>
	 * @throws NullPointerException
	 *             if {@code data} is null or {@code coeffs} is null
	 */
	public float[] smooth(float[] data, int from, int to) {
		return smooth(data, from, to, 0);
	}

	/**
	 * Runs filter on data from {@code from} (including) to {@code to}
	 * (excluding). Data beyond range spanned by {@code from} and {@code to}
	 * will be used for padding. See
	 * {@link #smooth(float[], float[], float[], int, float[][])} for usage
	 * of {@code bias}
	 * 
	 * @param data
	 *            data for filter
	 * @param from
	 *            inedx of the first element of data
	 * @param to
	 *            index of the first element omitted
	 * @param bias
	 *            how many points of pad should be left out when smoothing
	 * @param coeffs
	 *            filter coefficients
	 * @return filtered data
	 * @throws ArrayIndexOutOfBoundsException
	 *             if <tt>to > data.length</tt> or when {@code coeffs} has less
	 *             than <tt>2*bias + 1</tt> elements
	 * @throws IllegalArgumentException
	 *             if <tt>from < 0</tt> or <tt>to > data.length</tt> or
	 *             <tt>from > to</tt> or when <tt>bias < 0</tt> or
	 *             <tt>bias > min(nr, nl)</tt>
	 * @throws NullPointerException
	 *             if {@code data} is null or {@code coeffs} is null
	 */
	public float[] smooth(float[] data, int from, int to, int bias) {
		float[] leftPad=Arrays.copyOfRange(data, 0, from);
		float[] rightPad=Arrays.copyOfRange(data, to, data.length);
		float[] dataCopy=Arrays.copyOfRange(data, from, to);
		return smooth(dataCopy, leftPad, rightPad, bias);
	}
}
