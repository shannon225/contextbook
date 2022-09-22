package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;

import java.io.File;
import java.util.List;
import java.util.Optional;

public interface RetentionTimeAlignmentInterface extends ScoredPSMFilterInterface {
	float getYValue(float xrt);

	float getXValue(float yrt);

	float getProbabilityFitsModel(float actualRT, float modelRT);
	
	float getDelta(float actualRT, float modelRT);

	List<AlignmentDataPoint> plot(List<XYPoint> rts, Optional<File> saveFileSeed);

	interface AlignmentDataPoint {
		/**
		 * @return X value from library, in minutes
		 */
		float getLibrary();

		/**
		 * @return Y value from sample, in minutes
		 */
		float getActual();

		/**
		 * @return predicted Y value (in sample), in minutes
		 */
		float getPredictedActual();

		/**
		 * @return the smaller absolute delta (in minutes) from
		 *         this (x, y) point to the warping function
		 *         along the X or Y axes
		 */
		float getDelta();

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

		static AlignmentDataPoint of(float lib, float actual, float pred, float delta, float prob, Boolean decoy, String peptideModSeq) {
			return new AlignmentDataPoint() {
				@Override
				public float getLibrary() {
					return lib;
				}

				@Override
				public float getActual() {
					return actual;
				}

				@Override
				public float getPredictedActual() {
					return pred;
				}

				@Override
				public float getDelta() {
					return delta;
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

				@Override
				public String toString() {
					return String.format("(%.02f, %.02f, %.02f, %.02f, %s)",
							getLibrary(),
							getActual(),
							getPredictedActual(),
							getProbability(),
							getPeptideModSeq()
					);
				}
			};
		}
	}
}