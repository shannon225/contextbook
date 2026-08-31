package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

// Keeps only the scans acquired for the isolation window being scored
public class IsolationWindowFilter {
	private static final double EPSILON=0.01;
	private static final String PROPERTY="encyclopedia.exactIsolationWindow";

	public static boolean isEnabled() {
		return !"false".equalsIgnoreCase(System.getProperty(PROPERTY, "true"));
	}

	public static ArrayList<FragmentScan> restrictToWindow(ArrayList<FragmentScan> stripes, Range window) {
		if (!isEnabled()||stripes==null||stripes.isEmpty()||window==null) return stripes;

		if (countDistinctWindows(stripes)<=1) return stripes;

		ArrayList<FragmentScan> kept=new ArrayList<FragmentScan>(stripes.size());
		for (FragmentScan scan : stripes) {
			if (isFromWindow(scan, window)) {
				kept.add(scan);
			}
		}

		if (kept.isEmpty()) {
			Logger.logLine("No scan exactly matches "+window+" m/z; scoring against all "+stripes.size()
					+" overlapping scans.");
			return stripes;
		}

		Logger.logLine("Scheduled acquisition: "+window+" m/z is covered by scans from more than one isolation "
				+"window; keeping the "+kept.size()+" of "+stripes.size()+" acquired for this one. Disable with -D"
				+PROPERTY+"=false.");
		return kept;
	}

	private static int countDistinctWindows(ArrayList<FragmentScan> stripes) {
		Set<String> distinct=new HashSet<String>();
		for (FragmentScan scan : stripes) {
			distinct.add(windowKey(scan.getIsolationWindowLower(), scan.getIsolationWindowUpper()));
			if (distinct.size()>1) return distinct.size();
		}
		return distinct.size();
	}

	private static boolean isFromWindow(FragmentScan scan, Range window) {
		return Math.abs(scan.getIsolationWindowLower()-window.getStart())<EPSILON
				&&Math.abs(scan.getIsolationWindowUpper()-window.getStop())<EPSILON;
	}

	private static String windowKey(double lower, double upper) {
		return Math.round(lower/EPSILON)+"_"+Math.round(upper/EPSILON);
	}
}
