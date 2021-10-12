package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;

public class MultipleEnzymePrositCSVWriterTest {
	public static void main(String[] args) throws Exception {
		File fasta = new File("/Users/searleb/Documents/swaney/splice_junction_peptides/uniprot_human_reviewed_2021Oct6.fasta");
		File fastaWithSplicing = new File("/Users/searleb/Documents/swaney/splice_junction_peptides/uniprot_human_reviewed_plus_isoforms_2021Oct6.fasta");
		File justSplicing = new File("/Users/searleb/Documents/swaney/splice_junction_peptides/uniprot_human_reviewed_only_isoform_specific_2021Oct6.fasta");
		byte defaultNCE=30;
		byte defaultCharge=3;
		
		DigestionEnzyme[] enzymes=new DigestionEnzyme[] {DigestionEnzyme.getEnzyme("Trypsin"), DigestionEnzyme.getEnzyme("Glu-C"), DigestionEnzyme.getEnzyme("Asp-N")};
		for (int i = 0; i < enzymes.length; i++) {
			String fileName = PrositCSVWriter.checkCSVName(null, justSplicing, enzymes[i], defaultNCE, defaultCharge);
			
			HashSet<PeptidePrecursor> normalPeptides=PrositCSVWriter.getPeptidesFromFASTA(fasta, enzymes[i], (byte)2, (byte)3, 1, new Range(400, 1000));
			HashSet<PeptidePrecursor> splicingPeptides=PrositCSVWriter.getPeptidesFromFASTA(fastaWithSplicing, enzymes[i], (byte)2, (byte)3, 1, new Range(400, 1000));
			Logger.logLine("BEFORE: normalPeptides:"+normalPeptides.size()+", splicingPeptides:"+splicingPeptides.size());
			splicingPeptides.removeAll(normalPeptides);
			Logger.logLine("AFTER: normalPeptides:"+normalPeptides.size()+", splicingPeptides:"+splicingPeptides.size());
			

			int total = PrositCSVWriter.writePrositFile(fileName, defaultNCE, defaultCharge, false, splicingPeptides);
			Logger.logLine("Finished writing "+total+" peptides to Prosit CSV!");
		}
		
	}

}
