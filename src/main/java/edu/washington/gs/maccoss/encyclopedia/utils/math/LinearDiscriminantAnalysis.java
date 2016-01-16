package edu.washington.gs.maccoss.encyclopedia.utils.math;

public class LinearDiscriminantAnalysis implements ScoreCombiner {
	private final double[] coefficients;
	private final double constant;
	

	public LinearDiscriminantAnalysis(double[] coefficients, double constant) {
		this.coefficients=coefficients;
		this.constant=constant;
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

	public static LinearDiscriminantAnalysis buildModel(double[][] positiveData, double[][] negativeData) {
		double posPrior=positiveData.length/(float)(positiveData.length+negativeData.length);
		double negPrior=negativeData.length/(float)(positiveData.length+negativeData.length);

		int featureCount=Math.min(positiveData[0].length, negativeData[0].length);
		boolean[] useFeature=new boolean[featureCount];
		for (int i=0; i<useFeature.length; i++) {
			useFeature[i]=MatrixMath.getRange(MatrixMath.getColumn(positiveData, i))>0.0001&&MatrixMath.getRange(MatrixMath.getColumn(negativeData, i))>0.0001;
		}
		double[][] posData=selectUsedFeatures(positiveData, useFeature);
		double[][] negData=selectUsedFeatures(negativeData, useFeature);
		
		double[] meanAll=new double[featureCount];
		double[] meanPos=new double[featureCount];
		double[] meanNeg=new double[featureCount];
		double[] meanSum=new double[featureCount];
		double[] meanDiff=new double[featureCount];
		for (int i=0; i<meanAll.length; i++) {
			meanPos[i]=General.mean(MatrixMath.getColumn(posData, i));
			meanNeg[i]=General.mean(MatrixMath.getColumn(negData, i));
			meanAll[i]=meanPos[i]*posPrior+meanNeg[i]*negPrior;
			meanSum[i]=meanPos[i]+meanNeg[i];
			meanDiff[i]=meanPos[i]-meanNeg[i];
		}
		
		double[][] meanCorrectedPos=MatrixMath.subtract(posData, meanAll);
		double[][] meanCorrectedNeg=MatrixMath.subtract(negData, meanAll);
		
		double[][] covarPos=MatrixMath.multiply(MatrixMath.multiply(MatrixMath.transpose(meanCorrectedPos), meanCorrectedPos), 1.0/meanCorrectedPos.length);
		double[][] covarNeg=MatrixMath.multiply(MatrixMath.multiply(MatrixMath.transpose(meanCorrectedNeg), meanCorrectedNeg), 1.0/meanCorrectedNeg.length);

		double[][] pooledCovar=new double[covarPos.length][];
		for (int i=0; i<pooledCovar.length; i++) {
			pooledCovar[i]=new double[covarPos[i].length];
			for (int j=0; j<pooledCovar[i].length; j++) {
				pooledCovar[i][j]=covarPos[i][j]*posPrior+covarNeg[i][j]*negPrior;
			}
		}
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

