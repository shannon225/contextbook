package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import com.google.common.collect.ImmutableMap;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.impl.unmodifiable.TUnmodifiableObjectFloatMap;
import gnu.trove.map.TObjectFloatMap;
import gnu.trove.map.hash.TObjectFloatHashMap;

import java.util.Map;

/**
 * Interface for various implementations of intensity normalization.
 *
 * @see LibraryReportExtractor
 */
public interface IntensityNormalizer {
	float normalize(float intensity, String sourceFile);

	static IntensityNormalizer tic(TObjectFloatMap<? super String> ticBySourceFileMap) {
		// Immutable copy
		final TObjectFloatMap<? super String> ticMap = new TUnmodifiableObjectFloatMap<>(new TObjectFloatHashMap<>(ticBySourceFileMap));

		final float averageTic;
		if (ticMap.size() < 1) {
			averageTic = 0f;
		} else {
			averageTic = General.sum(ticMap.values()) / ticMap.size();
		}

		return (intensity, sourceFile) -> {
			float tic = ticMap.get(sourceFile);
			if (tic > 0.0f) {
				return intensity / tic * averageTic;
			} else {
				return intensity;
			}
		};
	}

	/**
	 * No normalization. Usable as {@code IntensityNormalizer::identity}.
	 */
	static float identity(float intensity, String ignored) {
		return intensity;
	}

	/**
	 * No normalization. Usable as {@code IntensityNormalizer.IDENTITY}.
	 */
	IntensityNormalizer IDENTITY = IntensityNormalizer::identity;

	/**
	 * No normalization. Usable as {@code IntensityNormalizer.identity()}.
	 */
	static IntensityNormalizer identity() {
		return IDENTITY;
	}
}
