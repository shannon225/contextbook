package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import uk.ac.ebi.pride.utilities.pridemod.ModReader;
import uk.ac.ebi.pride.utilities.pridemod.model.PTM;

public class ExtendedFastaEntry extends FastaEntry {
	private static ModReader modReader = ModReader.getInstance();
	private static HashMap<String, PTM> modMap=new HashMap<>();
	
	// immutable after construction
	private final ArrayList<AlleleVariant> potentialVariants=new ArrayList<AlleleVariant>();

	public ExtendedFastaEntry(String filename, String annotation, String sequence, SearchParameters parameters) {
		super(filename, annotation, sequence);
		String processedAnnotation=this.getAnnotation();
		String[] annotates=processedAnnotation.split(" \\\\");

		for (int index=0; index<annotates.length; index++) {
			if (parameters.getAAConstants().getVariableMods().isEmpty()) {
				// search for variants (don't do both)
				if (annotates[index].startsWith("VariantSimple=")) {
					parseVariantAnnotation(annotates[index].substring(13), true);
				} else if (annotates[index].startsWith("VariantComplex=")) {
					parseVariantAnnotation(annotates[index].substring(14), false);
				}
			} else {
				// otherwise search for mods (don't do both)
				if (annotates[index].startsWith("ModResPsi=")) {
					parseModificationAnnotation(annotates[index].substring(10), parameters);
				}
			}
		}
	}
	
	private void parseModificationAnnotation(String variantAnnotation, SearchParameters parameters) {
		String sequence=this.getSequence();
		int endIndex=0;
		for (int index=variantAnnotation.indexOf('('); index>=0; index=variantAnnotation.indexOf('(', endIndex+1)) {
			endIndex=variantAnnotation.indexOf(')', index);
			
			String thisMod = variantAnnotation.substring(index+1, endIndex);
			
			String[] info=thisMod.split("\\|");
			int start;
			try {
				start=Integer.parseInt(info[0]);
			} catch (NumberFormatException e) {
				throw new EncyclopediaException("Error on parsing variant index"+this.getAccession()+" ["+thisMod+"]", e);
			}

			char switchFrom=sequence.charAt(start-1);
			String modAccession=info[1];

			PTM ptm=modMap.get(modAccession);
			if (ptm==null) {
				ptm=modReader.getPTMbyAccession(modAccession);
				modMap.put(modAccession, ptm);
			}

			if (null==ptm) {
				Logger.errorLine("Unexpected modification ["+modAccession+"]!");
			}
			
			if (checkMod(parameters, switchFrom, ptm, start==1, start==sequence.length())) {
				potentialVariants.add(new AlleleVariant(start, ptm.getMonoDeltaMass(), switchFrom));
			}
		}
	}

	private boolean checkMod(SearchParameters parameters, char switchFrom, PTM ptm, boolean isN, boolean isC) {
		MassTolerance tolerance = parameters.getFragmentTolerance();
		ModificationMassMap variableMods = parameters.getAAConstants().getVariableMods();
		
		if (isN) {
			double delta = variableMods.getProteinNTermMod(switchFrom);
			if (tolerance.equals(ptm.getMonoDeltaMass(), delta)) return true;
		}
		if (isC) {
			double delta = variableMods.getProteinCTermMod(switchFrom);
			if (tolerance.equals(ptm.getMonoDeltaMass(), delta)) return true;
		}
		double delta = variableMods.getVariableMod(switchFrom);
		return tolerance.equals(ptm.getMonoDeltaMass(), delta);
	}

	/**
	 * only called during construction
	 * @param variantAnnotation
	 * @param simple
	 */
	private void parseVariantAnnotation(String variantAnnotation, Boolean simple) {
		String sequence=this.getSequence();
		int endIndex=0;
		for (int index=variantAnnotation.indexOf('('); index>=0; index=variantAnnotation.indexOf('(', endIndex+1)) {
			endIndex=variantAnnotation.indexOf(')', index);
			String[] info=variantAnnotation.substring(index+1, endIndex).split("\\|");
			int start;
			try {
				start=Integer.parseInt(info[0]);
			} catch (NumberFormatException e) {
				throw new EncyclopediaException("Error on parsing variant index"+this.getAccession()+" ["+variantAnnotation.substring(index+1, endIndex)+"]", e);
			}
			int end;
			if (simple) {
				char switchFrom=sequence.charAt(start-1);
				char switchTo=info[1].charAt(0);
				if (!(switchFrom=='I'&&switchTo=='L')&&!(switchFrom=='L'&&switchTo=='I')) {
					// if both aren't I/L then add
					potentialVariants.add(new AlleleVariant(start, switchFrom, switchTo));
				}
			} else {
				try {
					end=Integer.parseInt(info[1]);
				} catch (NumberFormatException e) {
					throw new EncyclopediaException("Error on parsing variant index"+this.getAccession()+" ["+variantAnnotation.substring(index+1, endIndex)+"]", e);
				}
				// length == 2 for the deletion cases, e.g. (4|4|)
				String newseq=(info.length>2)?info[2]:"";
				potentialVariants.add(new AlleleVariant(start, end, sequence.substring(start-1, end), newseq));
			}
		}
	}

	public ArrayList<AlleleVariant> getPotentialVariants() {
		return potentialVariants;
	}

}
