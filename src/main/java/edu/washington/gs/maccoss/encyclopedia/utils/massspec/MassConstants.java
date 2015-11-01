package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TIntCharHashMap;
import gnu.trove.procedure.TCharDoubleProcedure;

public class MassConstants {
	final private static TCharDoubleHashMap massesByAA=new TCharDoubleHashMap();
	final private static TIntCharHashMap aasByNominal=new TIntCharHashMap();
	private static void checkInit() {
		if (massesByAA.size()>0) return;
		
		massesByAA.put('A', 71.0371);
		massesByAA.put('R', 156.1011);
		massesByAA.put('N', 114.0429);
		massesByAA.put('D', 115.027);
		massesByAA.put('C', 103.0092);
		massesByAA.put('E', 129.0426);
		massesByAA.put('Q', 128.0586);
		massesByAA.put('G', 57.0215);
		massesByAA.put('H', 137.0589);
		massesByAA.put('L', 113.0841);
		massesByAA.put('I', 113.0841);
		massesByAA.put('K', 128.095);
		massesByAA.put('M', 131.0405);
		massesByAA.put('F', 147.0684);
		massesByAA.put('P', 97.0528);
		massesByAA.put('S', 87.032);
		massesByAA.put('T', 101.0477);
		massesByAA.put('W', 186.0793);
		massesByAA.put('Y', 163.0633);
		massesByAA.put('V', 99.0684);
		
		massesByAA.forEachEntry(new TCharDoubleProcedure() {
			public boolean execute(char arg0, double arg1) {
				aasByNominal.put((int)Math.round(arg1), arg0);
				return true;
			}
		});
	}
	
	public static double getMass(char aa) {
		checkInit();
		return massesByAA.get(aa);
	}
	
	public static double getMass(String sequence) {
		double total=0.0;
		for (char c : sequence.toCharArray()) {
			total+=getMass(c);
		}
		return total;
	}
	
	public static Float getModificationMass(String mod) {
		if ("Cam".equals(mod)) {
			return 57.0214635f;
		} else if ("O".equals(mod)) {
			return 15.994915f;
		}
		return null;
	}
	
	private static final MassTolerance tolerance=new MassTolerance(1.0); // 1 ppm is about the accuracy of floats 
	public static double getNeutralLoss(double modificationMass) {
		if (tolerance.equals(80.0, modificationMass)) {
			return 97.976896;
		} else if (tolerance.equals(79.966331, modificationMass)) {
			return 97.976896;
		}
		return 0.0;
	}
	
	public static Character getNearestAA(double mass) {
		checkInit();
		char c=aasByNominal.get((int)Math.round(mass));
		if (c!=0) {
			return c;
		} else {
			return null;
		}
	}
}
