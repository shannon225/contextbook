package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.procedure.TCharDoubleProcedure;

public class ModificationMassMap {
	private static final String NO_MODIFICATIONS="-";
	private static final String MOD_DELIMINATOR=",";
	private static final String EQUATION_DELIMINATOR="=";
	
	public static final char PROTEIN_C_TERM='z';
	public static final char PROTEIN_N_TERM='a';
	public static final char C_TERM='c';
	public static final char N_TERM='n';

	public static final double MISSING=0.0;

	private final boolean isEmpty;
	private final TCharDoubleHashMap primaryMassMap=new TCharDoubleHashMap();
	private final TCharDoubleHashMap ntermToMassMap=new TCharDoubleHashMap();
	private final TCharDoubleHashMap ctermToMassMap=new TCharDoubleHashMap();
	private final TCharDoubleHashMap proNtermToMassMap=new TCharDoubleHashMap();
	private final TCharDoubleHashMap proCtermToMassMap=new TCharDoubleHashMap();
	
	public ModificationMassMap() {
		// add no modifications
		isEmpty=true;
	}

	public ModificationMassMap(String value) {
		if (value!=null) {
			if (!NO_MODIFICATIONS.equals(value)) {
				try {
					StringTokenizer st=new StringTokenizer(value, MOD_DELIMINATOR);
					while (st.hasMoreTokens()) {
						String token=st.nextToken();
						StringTokenizer st2=new StringTokenizer(token, EQUATION_DELIMINATOR);
						String aa=st2.nextToken();
						double mass=Double.parseDouble(st2.nextToken());
						if (aa.length()==1) {
							primaryMassMap.put(aa.charAt(0), mass);
						} else if (aa.length()==2) {
							if (aa.charAt(0)==N_TERM) {
								ntermToMassMap.put(aa.charAt(1), mass);
							} else if (aa.charAt(0)==C_TERM) {
								ctermToMassMap.put(aa.charAt(1), mass);
							} else if (aa.charAt(0)==PROTEIN_N_TERM) {
								proNtermToMassMap.put(aa.charAt(1), mass);
							} else if (aa.charAt(0)==PROTEIN_C_TERM) {
								proCtermToMassMap.put(aa.charAt(1), mass);
							}
						} else {
							throw new EncyclopediaException("Error parsing variable modifications from ["+value+"], expected at most 2 characters, the terminal indicator and the amino acid");
						}
					}
				} catch (Exception e) {
					throw new EncyclopediaException("Error parsing variable modifications from ["+value+"]", e);
				}
			}
		}
		if (primaryMassMap.size()>0||ntermToMassMap.size()>0||ctermToMassMap.size()>0||proNtermToMassMap.size()>0||proCtermToMassMap.size()>0) {
			isEmpty=false;
		} else {
			isEmpty=true;
		}
	}
	
	public boolean isEmpty() {
		return isEmpty;
	}

	public String toString() {
		if (isEmpty==true) return NO_MODIFICATIONS;
		
		final ArrayList<String> mods=new ArrayList<String>();
		primaryMassMap.forEachEntry(new TCharDoubleProcedure() {
			@Override
			public boolean execute(char arg0, double arg1) {
				StringBuilder sb=new StringBuilder();
				sb.append(arg0);
				sb.append(EQUATION_DELIMINATOR);
				sb.append(arg1);
				mods.add(sb.toString());
				return true;
			}
		});
		ntermToMassMap.forEachEntry(new TCharDoubleProcedure() {
			@Override
			public boolean execute(char arg0, double arg1) {
				StringBuilder sb=new StringBuilder();
				sb.append(N_TERM);
				sb.append(arg0);
				sb.append(EQUATION_DELIMINATOR);
				sb.append(arg1);
				mods.add(sb.toString());
				return true;
			}
		});
		ctermToMassMap.forEachEntry(new TCharDoubleProcedure() {
			@Override
			public boolean execute(char arg0, double arg1) {
				StringBuilder sb=new StringBuilder();
				sb.append(C_TERM);
				sb.append(arg0);
				sb.append(EQUATION_DELIMINATOR);
				sb.append(arg1);
				mods.add(sb.toString());
				return true;
			}
		});
		proNtermToMassMap.forEachEntry(new TCharDoubleProcedure() {
			@Override
			public boolean execute(char arg0, double arg1) {
				StringBuilder sb=new StringBuilder();
				sb.append(PROTEIN_N_TERM);
				sb.append(arg0);
				sb.append(EQUATION_DELIMINATOR);
				sb.append(arg1);
				mods.add(sb.toString());
				return true;
			}
		});
		proCtermToMassMap.forEachEntry(new TCharDoubleProcedure() {
			@Override
			public boolean execute(char arg0, double arg1) {
				StringBuilder sb=new StringBuilder();
				sb.append(PROTEIN_C_TERM);
				sb.append(arg0);
				sb.append(EQUATION_DELIMINATOR);
				sb.append(arg1);
				mods.add(sb.toString());
				return true;
			}
		});
		Collections.sort(mods);

		StringBuilder sb=new StringBuilder();
		for (String mod : mods) {
			if (sb.length()>0) {
				sb.append(MOD_DELIMINATOR);
			}
			sb.append(mod);
		}
		return sb.toString();
	}

	public double getVariableMod(char c) {
		if (isEmpty==true) return MISSING;
		double value=primaryMassMap.get(c);
		return value;
	}

	public double getNTermMod(char c) {
		if (isEmpty==true) return MISSING;
		double value=primaryMassMap.get('n');
		if (value!=MISSING) return value;

		return ntermToMassMap.get(c);
	}

	public double getProteinNTermMod(char c) {
		if (isEmpty==true) return MISSING;
		double value=primaryMassMap.get('a');
		if (value!=MISSING) return value;

		return proNtermToMassMap.get(c);
	}

	public double getCTermMod(char c) {
		if (isEmpty==true) return MISSING;
		double value=primaryMassMap.get('c');
		if (value!=MISSING) return value;

		return ctermToMassMap.get(c);
	}

	public double getProteinCTermMod(char c) {
		if (isEmpty==true) return MISSING;
		double value=primaryMassMap.get('z');
		if (value!=MISSING) return value;

		return proCtermToMassMap.get(c);
	}
}
