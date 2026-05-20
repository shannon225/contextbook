package edu.washington.gs.maccoss.encyclopedia.context;

public class IsolationWindow {
	private String compound;
	private double targetMz;
	private byte charge;
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

	// Constructor
	public IsolationWindow(double precursorMz, byte charge, float rtInSecondsStart, float rtInSecondsStop, boolean isDecoy) {
		this.targetMz = precursorMz;
		this.charge = charge;
		this.rtMin = rtInSecondsStart;
		this.rtMax = rtInSecondsStop;
		this.isDecoy = isDecoy;
	}

	// Constructor
	public IsolationWindow(String compound, double precursorMz, byte charge, float rtInSecondsStart, float rtInSecondsStop, boolean isDecoy) {
		this.compound = compound;
		this.targetMz = precursorMz;
		this.charge = charge;
		this.rtMin = rtInSecondsStart;
		this.rtMax = rtInSecondsStop;
		this.isDecoy = isDecoy;
	}


	// Getters 
	
	public String getCompound() {
		return compound;
	}
	public double getTargetMz() {
		return targetMz;
	}

	public byte getCharge() {
		return charge;
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

