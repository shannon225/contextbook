package edu.washington.gs.maccoss.encyclopedia.context;

public class IsolationWindow {
	private double targetMz;
	private double mzStart;
	private double mzStop; 
	private double windowMz;
	private float rtCenter; // Including this because we use it to calculate m/z min and max in the TargetedBoostrapper
	private float rtMin;
	private float rtMax;
	private boolean isDecoy;

	
	// Constructor
	public IsolationWindow(double precursorMz, float rtInSecondsStart, float rtInSecondsStop, boolean isDecoy) {
		this.targetMz = precursorMz;
		this.rtMin = rtInSecondsStart;
		this.rtMax = rtInSecondsStop;
		this.isDecoy = isDecoy;
		}


	// Getters 
	public double getTargetMz() {
		return targetMz;
	}
	
	public double getWindowMz() {
		return windowMz;
	}
	
	public double getMzStart() {
		return mzStart;
	}
	
	public double getMzStop() {
		return mzStop;
	}

	public float getRtCenter() {
		return rtCenter;
	}
	
	public float getRtMin() {
		return rtMin;
	}
	
	public float getRtMax() {
		return rtMax;
	}
	
	public boolean isDecoy() {
		return isDecoy;
	}

	public String size() {
		// TODO Auto-generated method stub
		return null;
	}

}

