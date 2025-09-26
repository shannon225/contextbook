package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class LinearDiscriminantAnalysis implements ScoreCombiner, Comparable<LinearDiscriminantAnalysis> {
	private final double[] coefficients;
	private final double constant;
	

	public LinearDiscriminantAnalysis(double[] coefficients, double constant) {
		this.coefficients=coefficients;
		this.constant=constant;
	}
	
	@Override
	public int compareTo(LinearDiscriminantAnalysis o) {
		if (o==null) return 1;
		// LDAs aren't directly comparable, so this is just to maintain consistent sort order!
		return Double.compare(General.sum(coefficients), General.sum(o.coefficients));
	}
	
	public static LinearDiscriminantAnalysis average(ArrayList<LinearDiscriminantAnalysis> ldas) {
		double[] coefficients=null;
		double constant=0.0;
		
		int count=0;
		for (LinearDiscriminantAnalysis lda : ldas) {
			if (coefficients==null) {
				coefficients=lda.coefficients.clone();
				constant=lda.constant;
			} else {
				coefficients=MatrixMath.add(coefficients, lda.coefficients);
				constant+=lda.constant;
			}
			count++;
		}
		
		coefficients=General.divide(coefficients, count);
		constant=constant/count;
		return new LinearDiscriminantAnalysis(coefficients, constant);
	}
	
	@Override
	public String toString() {
		StringBuilder sb=new StringBuilder("coeff:[");
		for (int i=0; i<coefficients.length; i++) {
			if (i>0) sb.append(", ");
			sb.append(coefficients[i]);
		}
		sb.append("], c:");
		sb.append(constant);
		return sb.toString();
	}

	@Override
	public float getScore(float[] data) {
		int length=Math.min(data.length, coefficients.length);
		double score=constant;
		for (int i=0; i<length; i++) {
			score+=data[i]*coefficients[i];
		}
		return (float)score;
	}
	
	public double[] getCoefficients() {
		return coefficients;
	}
	
	public double getConstant() {
		return constant;
	}
	
	public static LinearDiscriminantAnalysis buildModel(ArrayList<float[]> positiveData, ArrayList<float[]> negativeData) {
		double[][] positive=new double[positiveData.size()][];
		double[][] negative=new double[negativeData.size()][];
		
		for (int i=0; i<positive.length; i++) {
			positive[i]=General.toDoubleArray(positiveData.get(i));
		}
		for (int i=0; i<negative.length; i++) {
			negative[i]=General.toDoubleArray(negativeData.get(i));
		}
		return buildModel(positive, negative);
	}

	public static LinearDiscriminantAnalysis buildModel(double[][] positiveData, double[][] negativeData) {
		try {
		double posPrior=positiveData.length/(float)(positiveData.length+negativeData.length);
		double negPrior=negativeData.length/(float)(positiveData.length+negativeData.length);

		int featureCount=Math.min(positiveData[0].length, negativeData[0].length);
		boolean[] useFeature=new boolean[featureCount];
		for (int i=0; i<useFeature.length; i++) {
			useFeature[i]=MatrixMath.getRange(MatrixMath.getColumn(positiveData, i))>0.0001&&MatrixMath.getRange(MatrixMath.getColumn(negativeData, i))>0.0001;
		}
		double[][] posData=selectUsedFeatures(positiveData, useFeature);
		double[][] negData=selectUsedFeatures(negativeData, useFeature);
		featureCount=Math.min(posData[0].length, negData[0].length);
		
		double[] meanAll=new double[featureCount];
		double[] meanPos=new double[featureCount];
		double[] meanNeg=new double[featureCount];
		double[] meanSum=new double[featureCount];
		double[] meanDiff=new double[featureCount];
		for (int i=0; i<meanAll.length; i++) {
			meanPos[i]=General.mean(MatrixMath.getColumn(posData, i));
			meanNeg[i]=General.mean(MatrixMath.getColumn(negData, i));
			
			// meanAll uses the priors to calculate the weighted mean
			meanAll[i]=meanPos[i]*posPrior+meanNeg[i]*negPrior;
			
			meanSum[i]=meanPos[i]+meanNeg[i];
			meanDiff[i]=meanPos[i]-meanNeg[i];
		}
		
		// normalize by the average of the columns in both matrices 
		double[][] meanCorrectedPos=MatrixMath.subtract(posData, meanAll);
		double[][] meanCorrectedNeg=MatrixMath.subtract(negData, meanAll);
		
		// calculate the covariance matrix for both positive and negative matrices
		double[][] covarPos=MatrixMath.multiply(MatrixMath.multiply(MatrixMath.transpose(meanCorrectedPos), meanCorrectedPos), 1.0/meanCorrectedPos.length);
		double[][] covarNeg=MatrixMath.multiply(MatrixMath.multiply(MatrixMath.transpose(meanCorrectedNeg), meanCorrectedNeg), 1.0/meanCorrectedNeg.length);

		// calculate the average covariance matrix, again using priors to weight
		double[][] pooledCovar=new double[covarPos.length][];
		for (int i=0; i<pooledCovar.length; i++) {
			pooledCovar[i]=new double[covarPos[i].length];
			for (int j=0; j<pooledCovar[i].length; j++) {
				pooledCovar[i][j]=covarPos[i][j]*posPrior+covarNeg[i][j]*negPrior;
				if (!Double.isFinite(pooledCovar[i][j])) {
					System.out.println(pooledCovar[i][j]+" = "+covarPos[i][j]+"*"+posPrior+" + "+covarNeg[i][j]+"*"+negPrior);
				}
			}
		}
		
		// The inverse covariance is a measure of precision, while covariance is a measure of dispersion. 
		// Increased dispersion occurs when values are farther apart and the more they co-vary with other 
		// variables. LDA wants to upweight variables that have low variance and do not co-vary with others. 
		double[][] inversePooledCovar=MatrixMath.invert(pooledCovar);
		
		double[] coefficients=MatrixMath.multiply(inversePooledCovar, meanDiff);
		
		double zeroPoint=MatrixMath.multiply(coefficients, meanSum);
		double constant=-Math.log(negPrior/posPrior)-0.5*zeroPoint;
		
		double[] allCoefficients=new double[useFeature.length];
		int index=0;
		for (int i=0; i<useFeature.length; i++) {
			if (!useFeature[i]) {
				allCoefficients[i]=0.0;
			} else {
				allCoefficients[i]=coefficients[index];
				index++;
			}
		}
		return new LinearDiscriminantAnalysis(allCoefficients, constant);
		} catch (ArrayIndexOutOfBoundsException aiobe) {
			Logger.errorLine("LDA indexing error. Starting data was pos: "+positiveData.length+" neg: "+negativeData.length);
			for (int i = 1; i < negativeData.length; i++) {
				if (negativeData[i].length!=negativeData[i-1].length) {
					Logger.errorLine("Inconsistent negative record length: "+negativeData[i].length+" != "+negativeData[i-1].length);
				}
			}
			Logger.errorLine("Reference negative record length: "+negativeData[0].length);
			for (int i = 1; i < positiveData.length; i++) {
				if (positiveData[i].length!=positiveData[i-1].length) {
					Logger.errorLine("Inconsistent positive record length: "+positiveData[i].length+" != "+positiveData[i-1].length);
				}
			}
			Logger.errorLine("Reference positive record length: "+positiveData[0].length);
			throw new EncyclopediaException("Error calculating LDA", aiobe);
		}

	}

	private static double[][] selectUsedFeatures(double[][] data, boolean[] useFeature) {
		double[][] negData=new double[data.length][];
		int usedFeatureCount=0;
		for (int i=0; i<useFeature.length; i++) {
			if (useFeature[i]) usedFeatureCount++;
		}
		for (int i=0; i<data.length; i++) {
			negData[i]=new double[usedFeatureCount];
			int index=0;
			for (int j=0; j<useFeature.length; j++) {
				if (useFeature[j]) {
					negData[i][index]=data[i][j];
					index++;
				}
			}
		}
		return negData;
	}

}

