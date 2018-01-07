package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class ExtendedFastaEntry extends FastaEntry {
	private final ArrayList<AlleleVariant> potentialVariants=new ArrayList<AlleleVariant>();

	public ExtendedFastaEntry(String filename, String annotation, String sequence) {
		super(filename, annotation, sequence);
		String processedAnnotation=this.getAnnotation();
		String[] annotates=processedAnnotation.split(" ");

		for (int index=0; index<annotates.length; index++) {
			if (annotates[index].startsWith("\\VariantSimple=")) {
				parseVariantAnnotation(annotates[index].substring(14), true);
			} else if (annotates[index].startsWith("\\VariantComplex=")) {
				parseVariantAnnotation(annotates[index].substring(15), false);
			}
		}
	}

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
				this.addPotentialVariant(new AlleleVariant(start, sequence.charAt(start-1), info[1].charAt(0)));
			} else {
				try {
					end=Integer.parseInt(info[1]);
				} catch (NumberFormatException e) {
					throw new EncyclopediaException("Error on parsing variant index"+this.getAccession()+" ["+variantAnnotation.substring(index+1, endIndex)+"]", e);
				}
				// length == 2 for the deletion cases, e.g. (4|4|)
				String newseq=(info.length>2)?info[2]:"";
				this.addPotentialVariant(new AlleleVariant(start, end, sequence.substring(start-1, end), newseq));
			}
		}
	}

	public void addPotentialVariant(AlleleVariant variant) {
		potentialVariants.add(variant);
	}

	public ArrayList<AlleleVariant> getPotentialVariant() {
		return potentialVariants;
	}

}
