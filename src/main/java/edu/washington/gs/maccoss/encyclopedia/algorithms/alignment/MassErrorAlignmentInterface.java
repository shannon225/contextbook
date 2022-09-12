package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;

public interface MassErrorAlignmentInterface {
	List<MassErrorDataPoint> plot(ArrayList<XYPoint> errorByRT, Optional<File> saveFileSeed);

	float getYValue(float xrt); // y is Error

	float getXValue(float yrt); // x is RT

	float getProbabilityFitsModel(float actualRT, float massError);
	
	float getCorrectedMassError(float actualRT, float massError);

	interface MassErrorDataPoint {
		/**
		 * @return X value from actualRT, in minutes
		 */
		float getRetentionTime();

		/**
		 * @return Y value from massError
		 */
		float getMassError();

		/**
		 * @return predicted Y value
		 */
		float getPredictedMassError();

		/**
		 * @return the absolute error
		 */
		float getCorrectedMassError();

		/**
		 * @return the post-hoc probability of this data point being
		 *         a true alignment given the computed alignment
		 */
		float getProbability();

		/**
		 * @return if the data point corresponds to decoy data
		 */
		Boolean isDecoy();

		/**
		 *
		 * @return the peptide mod sequence
		 */
		String getPeptideModSeq();

		static MassErrorDataPoint of(float retentionTime, float massError, float predictedMassError, float correctedAbsoluteMassError, float prob, Boolean decoy, String peptideModSeq) {
			return new MassErrorDataPoint() {
				

				@Override
				public float getRetentionTime() {
					return retentionTime;
				}

				@Override
				public float getMassError() {
					return massError;
				}

				@Override
				public float getPredictedMassError() {
					return predictedMassError;
				}

				@Override
				public float getCorrectedMassError() {
					return correctedAbsoluteMassError;
				}

				@Override
				public float getProbability() {
					return prob;
				}

				@Override
				public Boolean isDecoy() {
					return decoy;
				}

				@Override
				public String getPeptideModSeq() {
					return peptideModSeq;
				}
			};
		}
	}
}