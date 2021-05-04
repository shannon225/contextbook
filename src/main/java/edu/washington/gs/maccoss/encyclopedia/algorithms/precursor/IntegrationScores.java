package edu.washington.gs.maccoss.encyclopedia.algorithms.precursor;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PeakTrace;
import edu.washington.gs.maccoss.encyclopedia.datastructures.RetentionTimeBoundary;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PrecursorIon;

public class IntegrationScores {
	private final PrecursorIon ion;
	private final float intensity;
	private final float apexRetentionTime;
	private final float averagePeakCorrelation;
	private final float isotopicRatio;
	private final RetentionTimeBoundary boundary;
	private final PeakTrace<String> medianPeakShape;

	public IntegrationScores(PrecursorIon ion, float intensity, float apexRetentionTime, float averagePeakCorrelation,
			float isotopicRatio, RetentionTimeBoundary boundary, PeakTrace<String> medianPeakShape) {
		this.ion = ion;
		this.intensity = intensity;
		this.apexRetentionTime = apexRetentionTime;
		this.averagePeakCorrelation = averagePeakCorrelation;
		this.isotopicRatio = isotopicRatio;
		this.boundary=boundary;
		this.medianPeakShape=medianPeakShape;
	}

	public PrecursorIon getIon() {
		return ion;
	}

	public float getIntensity() {
		return intensity;
	}

	public float getApexRetentionTime() {
		return apexRetentionTime;
	}

	public float getAveragePeakCorrelation() {
		return averagePeakCorrelation;
	}

	public float getIsotopicRatio() {
		return isotopicRatio;
	}

	public RetentionTimeBoundary getBoundary() {
		return boundary;
	}
	
	public PeakTrace<String> getMedianPeakShape() {
		return medianPeakShape;
	}
}
