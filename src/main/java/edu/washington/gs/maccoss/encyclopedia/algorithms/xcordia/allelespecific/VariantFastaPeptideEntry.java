package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;

public class VariantFastaPeptideEntry extends FastaPeptideEntry {
	private final AlleleVariant variant;

	public VariantFastaPeptideEntry(FastaPeptideEntry peptide, AlleleVariant variant) {
		super(peptide.getFilename(), peptide.getAccessions(), peptide.getSequence());
		this.variant=variant;
	}
	public AlleleVariant getVariant() {
		return variant;
	}
}
