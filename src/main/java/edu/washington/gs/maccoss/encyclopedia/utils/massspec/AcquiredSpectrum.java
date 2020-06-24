package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

public interface AcquiredSpectrum extends Spectrum {
	public float getIonInjectionTime();
	
	public int getFraction();
	public double getIsolationWindowLower();
	public double getIsolationWindowUpper();
}
