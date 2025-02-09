package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashMap;

public class GlobalRangeTracker {
//	HashMap<Range, Range> precursorRTsInSecs=new HashMap<Range, Range>();
	HashMap<Range, Range> stripeRTsInSecs=new HashMap<Range, Range>();
	
	public HashMap<Range, Range> getStripeRTsInSecs() {
		return stripeRTsInSecs;
	}
	
//	public HashMap<Range, Range> getPrecursorRTsInSecs() {
//		return precursorRTsInSecs;
//	}
	
	/**
	 * NOT THREAD SAFE!
	 * @param mzRange
	 * @param rtRange
	 * @return returns if should continue
	 */
	public void addRange(Range mzRange, Range rtRange) {
		stripeRTsInSecs.put(mzRange, rtRange);
	}
	
//	/**
//	 * NOT THREAD SAFE!
//	 * @param mzRange
//	 * @param rtRange
//	 */
//	public void addPrecursor(Range mzRange, Range rtRange) {
//		precursorRTsInSecs.put(mzRange, rtRange);
//	}
}
