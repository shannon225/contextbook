package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TIntCharHashMap;
import gnu.trove.procedure.TCharDoubleProcedure;

public class MassConstants {
	public final static double neutronMass=1.0086649158849;
	public final static double protonMass=1.00727646681290;
	public final static double hydrogenMass=1.007825032071;
	public final static double oh3=15.9949146195616+3*hydrogenMass;
	
	final private static TCharDoubleHashMap massesByAA=new TCharDoubleHashMap();
	final private static TIntCharHashMap aasByNominal=new TIntCharHashMap();
	private static void checkInit() {
		if (massesByAA.size()>0) return;
		
		massesByAA.put('A', 71.037114);
		massesByAA.put('R', 156.101111);
		massesByAA.put('N', 114.042927);
		massesByAA.put('D', 115.026943);
		massesByAA.put('C', 103.009185);
		massesByAA.put('E', 129.042593);
		massesByAA.put('Q', 128.058578);
		massesByAA.put('G', 57.021464);
		massesByAA.put('H', 137.058912);
		massesByAA.put('L', 113.084064);
		massesByAA.put('I', 113.084064);
		massesByAA.put('K', 128.094963);
		massesByAA.put('M', 131.040485);
		massesByAA.put('F', 147.068414);
		massesByAA.put('P', 97.052764);
		massesByAA.put('S', 87.032028);
		massesByAA.put('T', 101.047679);
		massesByAA.put('W', 186.079313);
		massesByAA.put('Y', 163.06332);
		massesByAA.put('V', 99.068414);
		
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
	
	public static double getChargedMass(String sequence, byte charge) {
		double mass=getMass(sequence)+oh3;
		return (mass+protonMass*charge)/charge;
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
