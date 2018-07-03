package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.Collection;
import java.util.Comparator;
import java.util.Optional;

import com.google.common.collect.ImmutableCollection;
import com.google.common.collect.ImmutableList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TCharObjectHashMap;
import gnu.trove.map.hash.TIntCharHashMap;
import gnu.trove.procedure.TCharDoubleProcedure;

public class AminoAcidConstants {
	public static final char[] AAs="ARNDCEQGHLIKMFPSTWYV".toCharArray();
	
	// ordered by H C O N S
	private final TCharDoubleHashMap fixedMods;
	private final ModificationMassMap variableMods;
	private final ImmutableCollection<PeptideModification> localizationModifications;
	private final TCharObjectHashMap<int[]> atomicComposition=new TCharObjectHashMap<int[]>();
	final private TCharDoubleHashMap massesByAA=new TCharDoubleHashMap();
	final private TIntCharHashMap aasByNominal=new TIntCharHashMap();

	public static TCharDoubleHashMap getFixedModsMap(String name) {
		TCharDoubleHashMap fixedMods;
		if ("C+57 (Carbamidomethyl)".equalsIgnoreCase(name)) {
			fixedMods=new TCharDoubleHashMap(new char[] {'C'}, new double[] {57.0214635});
			
		} else if ("C+58 (Carboxymethyl)".equalsIgnoreCase(name)) {
			fixedMods=new TCharDoubleHashMap(new char[] {'C'}, new double[] {58.005479});
			
		}else if ("C+46 (MMTS)".equalsIgnoreCase(name)) {
			fixedMods=new TCharDoubleHashMap(new char[] {'C'}, new double[] {45.987721});
			
		} else {
			fixedMods=new TCharDoubleHashMap();
		}
		return fixedMods;
	}
	
	public static String toName(AminoAcidConstants constants) {
		if (Math.round(constants.getFixedMods().get('C'))==57) {
			return "C+57 (Carbamidomethyl)";
		} else if (Math.round(constants.getFixedMods().get('C'))==58) {
			return "C+58 (Carboxymethyl)";
		} else if (Math.round(constants.getFixedMods().get('C'))==46) {
			return "C+46 (MMTS)";
		}
		return "No fixed modifications";
	}

	public static ImmutableCollection<PeptideModification> getDefaultLocalizationModifications() {
		return ImmutableList.of(
				PeptideModification.phosphorylation,
				PeptideModification.acetylation,
				PeptideModification.oxidation,
				PeptideModification.methylation,
				PeptideModification.ubiquitination,
				PeptideModification.oHexNAc);
	}

	public static AminoAcidConstants getConstants(String fixedAAName, ModificationMassMap variableMods) {
		TCharDoubleHashMap fixedMods = getFixedModsMap(fixedAAName);
		return new AminoAcidConstants(fixedMods, variableMods);
	}

	public static AminoAcidConstants createEmptyFixedAndVariable() {
		return new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
	}

	/**
	 * assumes +57 C-alkylation
	 */
	public AminoAcidConstants() {
		this(new TCharDoubleHashMap(new char[] {'C'}, new double[] {57.0214635}), new ModificationMassMap());
	}

	public AminoAcidConstants(TCharDoubleHashMap fixedMods, ModificationMassMap variableMods) {
		this(fixedMods, variableMods, getDefaultLocalizationModifications());
	}

	public AminoAcidConstants(TCharDoubleHashMap fixedMods, ModificationMassMap variableMods, Collection<PeptideModification> localizationModifications) {
		this.fixedMods=fixedMods;
		this.variableMods=variableMods;
		this.localizationModifications = ImmutableList.copyOf(localizationModifications);

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
		
		for (int i=0; i<AAs.length; i++) {
			aasByNominal.put((int)Math.round(massesByAA.get(AAs[i])), AAs[i]);
		}
	}
	
	public TCharDoubleHashMap getFixedMods() {
		return fixedMods;
	}
	
	public ModificationMassMap getVariableMods() {
		return variableMods;
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
	
	public String getVariableModString() {
		return variableMods.toString();
	}
	
	public double getMass(char aa) {
		return massesByAA.get(aa);
	}
	
	public double getMass(String sequence) {
		double total=0.0;
		for (int i=0; i<sequence.length(); i++) {
			char c=sequence.charAt(i);

			if (c=='[') {
				final StringBuilder sb=new StringBuilder();
				while (true) {
					i++;
					c=sequence.charAt(i);
					if (c==']') {
						break; 
					} else {
						sb.append(c);
					}
				}
				final double modMass = Double.parseDouble(sb.toString());
				total += modMass;
			}
			total+=getMass(c);
		}
		return total;
	}
	
	public double getChargedMass(String modSeq, byte charge) {
		double mass=getMass(modSeq)+MassConstants.oh2;
		return (mass+MassConstants.protonMass*charge)/charge;
	}

	public double getChargedIsotopeMass(String modSeq, byte charge, byte isotope) {
		return MassConstants.getChargedIsotopeMass(getChargedMass(modSeq, charge), charge, isotope);
	}

	public int[] getAminoAcidProportions(char c) {
		int[] is=atomicComposition.get(c);
		return is;
	}
	
	public Character getNearestAA(double mass) {
		char c=aasByNominal.get((int)Math.round(mass));
		if (c!=0) {
			return c;
		} else {
			return null;
		}
	}

	public ImmutableCollection<PeptideModification> getLocalizationModifications() {
		return localizationModifications;
	}

	public double getNeutralLoss(char aminoAcid, double modificationMass) {
		return getNeutralLoss(localizationModifications, aminoAcid, modificationMass)
				.orElseGet(() -> getNeutralLoss(getDefaultLocalizationModifications(), aminoAcid, modificationMass).orElse(0d));
	}

	private static Optional<Double> getNeutralLoss(Collection<PeptideModification> modifications, char aminoAcid, double modificationMass) {
		return findModification(modifications, aminoAcid, modificationMass)
				.map(mod -> mod.getNeutralLoss(aminoAcid));
	}

	private static Optional<PeptideModification> findModification(Collection<PeptideModification> modifications, char aminoAcid, double modificationMass) {
		return modifications.stream()
				.filter(mod -> mod.isModificationMass(aminoAcid, modificationMass))
				.sorted(Comparator.comparing(mod -> Math.abs(mod.getMass() - modificationMass)))
				.findFirst();
	}

	private static final MassTolerance tolerance=new MassTolerance(1.0); // 1 ppm is about the accuracy of floats

	public double getAccurateModificationMass(char aa, double modificationMass) {

		Double fixedOrVariable;

		final double noEntry = fixedMods.getNoEntryValue();
		double fixedModMass = fixedMods.get(aa);
		if (noEntry == fixedModMass) {
			double variableModMass = variableMods.getVariableMod(aa);
			fixedOrVariable = variableModMass == ModificationMassMap.MISSING ? null : variableModMass;
		} else {
			fixedOrVariable = fixedModMass;
		}

		if (fixedOrVariable != null) {
			return fixedOrVariable;
		} else {
			Optional<PeptideModification> modification = findModification(getDefaultLocalizationModifications(), aa, modificationMass);
			Optional<Double> opModificationMass = modification.map(PeptideModification::getMass);
			return opModificationMass.orElseGet(() -> {
				if (aa == 'C') {
					if (tolerance.equals(57.0, modificationMass)) { // Carbamidomethyl
						return 57.0214635;
					} else if (tolerance.equals(58.0, modificationMass)) { // Carboxymethyl
						return 58.005479;
					} else if (tolerance.equals(46.0, modificationMass)) { // MMTS
						return 45.987721;
					} else if (tolerance.equals(99.0, modificationMass)) { // Carbamidomethyl + acetyl
						return 57.0214635 + 42.010565;
					} else if (tolerance.equals(40.0, modificationMass)) { // Carbamidomethyl - pyro-glu
						return 57.0214635 - 17.026549;
					}
				}

				if (aa == 'M' || aa == 'W') {
					if (tolerance.equals(16.0, modificationMass)) { // Oxidation
						return 15.994915;
					} else if (tolerance.equals(58.0, modificationMass)) { // Ox + acetyl
						return 42.010565 + 15.994915;
					}
				}

				if (aa == 'Q') {
					if (tolerance.equals(-17.0, modificationMass)) { // pyro-glu
						return -17.026549;
					}
				}

				if (aa == 'S' || aa == 'T' || aa == 'Y') {
					if (tolerance.equals(80.0, modificationMass)) { // Phospho
						return PeptideModification.phosphorylation.getMass();
					} else if (tolerance.equals(122.0, modificationMass)) { // Phospho + acetyl
						return 42.010565 + PeptideModification.phosphorylation.getMass();
					}
				}

				if (tolerance.equals(42.0, modificationMass)) { // acetyl
					return 42.010565;
				} else {
					return modificationMass;
				}
			});
		}
	}
}
