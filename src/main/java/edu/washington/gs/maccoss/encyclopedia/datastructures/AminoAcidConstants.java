package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TCharObjectHashMap;
import gnu.trove.map.hash.TIntCharHashMap;
import gnu.trove.procedure.TCharDoubleProcedure;

public class AminoAcidConstants {
	
	// ordered by H C O N S
	private final TCharDoubleHashMap fixedMods;
	private final TCharObjectHashMap<int[]> atomicComposition=new TCharObjectHashMap<int[]>();
	final private TCharDoubleHashMap massesByAA=new TCharDoubleHashMap();
	final private TIntCharHashMap aasByNominal=new TIntCharHashMap();
	
	public static AminoAcidConstants getConstants(String name) {
		if ("C+57 (Carbamidomethyl)".equalsIgnoreCase(name)) {
			return new AminoAcidConstants(new TCharDoubleHashMap(new char[] {'C'}, new double[] {57.0214635}));
			
		} else if ("C+58 (Carboxymethyl)".equalsIgnoreCase(name)) {
			return new AminoAcidConstants(new TCharDoubleHashMap(new char[] {'C'}, new double[] {58.005479}));
			
		}else if ("C+46 (MMTS)".equalsIgnoreCase(name)) {
			return new AminoAcidConstants(new TCharDoubleHashMap(new char[] {'C'}, new double[] {45.987721}));
			
		} else {
			return new AminoAcidConstants(new TCharDoubleHashMap());
		}
	}

	/**
	 * assumes +57 C-alkylation
	 */
	public AminoAcidConstants() {
		this(new TCharDoubleHashMap(new char[] {'C'}, new double[] {57.0214635}));
	}
	public AminoAcidConstants(TCharDoubleHashMap fixedMods) {
		this.fixedMods=fixedMods;
		
		atomicComposition.put('A', new int[] {5, 3, 1, 1, 0});
		if (fixedMods.contains('C')&&Math.round(fixedMods.get('C'))==57) {
			atomicComposition.put('C', new int[] {8, 5, 2, 2, 1}); // assumes +57 is carbamidomethyl alkylation
		} else if (fixedMods.contains('C')&&Math.round(fixedMods.get('C'))==58) {
			atomicComposition.put('C', new int[] {7, 5, 3, 1, 1}); // assumes +58 is carboxymethyl alkylation
		} else if (fixedMods.contains('C')&&Math.round(fixedMods.get('C'))==46) {
			atomicComposition.put('C', new int[] {7, 4, 1, 1, 2}); // assumes +46 is methylthio (MMTS)
		} else {
			atomicComposition.put('C', new int[] {5, 3, 1, 1, 1}); // unmodified
		}
		atomicComposition.put('D', new int[] {5, 4, 3, 1, 0});
		atomicComposition.put('E', new int[] {7, 5, 3, 1, 0});
		atomicComposition.put('F', new int[] {9, 9, 1, 1, 0});
		atomicComposition.put('G', new int[] {3, 2, 1, 1, 0});
		atomicComposition.put('H', new int[] {7, 6, 1, 3, 0});
		atomicComposition.put('I', new int[] {11, 6, 1, 1, 0});
		atomicComposition.put('K', new int[] {12, 6, 1, 2, 0});
		atomicComposition.put('L', new int[] {11, 6, 1, 1, 0});
		atomicComposition.put('M', new int[] {9, 5, 1, 1, 1});
		atomicComposition.put('N', new int[] {6, 4, 2, 2, 0});
		atomicComposition.put('P', new int[] {7, 5, 1, 1, 0});
		atomicComposition.put('Q', new int[] {8, 5, 2, 2, 0});
		atomicComposition.put('R', new int[] {12, 6, 1, 4, 0});
		atomicComposition.put('S', new int[] {5, 3, 2, 1, 0});
		atomicComposition.put('T', new int[] {7, 4, 2, 1, 0});
		atomicComposition.put('V', new int[] {9, 5, 1, 1, 0});
		atomicComposition.put('W', new int[] {10, 11, 1, 2, 0});
		atomicComposition.put('Y', new int[] {9, 9, 2, 1, 0});

		
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
		
		fixedMods.forEachEntry(new TCharDoubleProcedure() {
			@Override
			public boolean execute(char aa, double m) {
				massesByAA.adjustOrPutValue(aa, m, m);
				return true;
			}
		});
		
		massesByAA.forEachEntry(new TCharDoubleProcedure() {
			public boolean execute(char arg0, double arg1) {
				aasByNominal.put((int)Math.round(arg1), arg0);
				return true;
			}
		});
	}
	
	public TCharDoubleHashMap getFixedMods() {
		return fixedMods;
	}
	
	public String getFixedModString() {
		final StringBuilder sb=new StringBuilder();
		fixedMods.forEachEntry(new TCharDoubleProcedure() {
			@Override
			public boolean execute(char arg0, double arg1) {
				if (sb.length()>0) {
					sb.append(",");
				}
				sb.append(arg0);
				sb.append("=");
				sb.append(arg1);
				return true;
			}
		});
		return sb.toString();
	}
	
	public double getMass(char aa) {
		return massesByAA.get(aa);
	}
	
	public double getMass(String sequence) {
		double total=0.0;
		for (char c : sequence.toCharArray()) {
			total+=getMass(c);
		}
		return total;
	}
	
	public double getChargedMass(String sequence, byte charge) {
		double mass=getMass(sequence)+MassConstants.oh2;
		return (mass+MassConstants.protonMass*charge)/charge;
	}
	
	public int[] getAminoAcidProportions(char c) {
		return atomicComposition.get(c);
	}
	
	public Character getNearestAA(double mass) {
		char c=aasByNominal.get((int)Math.round(mass));
		if (c!=0) {
			return c;
		} else {
			return null;
		}
	}
}
