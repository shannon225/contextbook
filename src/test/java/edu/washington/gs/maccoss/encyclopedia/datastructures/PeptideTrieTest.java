package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.InputStream;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.filereaders.MSPReader;
import junit.framework.TestCase;

public class PeptideTrieTest extends TestCase {

	public void testTrie() throws Exception {
		InputStream is=getClass().getResourceAsStream("/truncated.msp");
		ArrayList<LibraryEntry> entries=MSPReader.readMSP(is, "truncated.msp", true);
		PeptideAccessionMatchingTrie trie=new PeptideAccessionMatchingTrie(entries) {
			@Override
			protected void processMatch(FastaEntryInterface fasta, LibraryEntry entry, int start) {
				super.processMatch(fasta, entry, start);
				PeptideTrieTest.assertEquals(entry.getPeptideSeq(), fasta.getSequence().substring(start, start+entry.getPeptideSeq().length()));
			}
		};
		
		FastaEntry fasta=new FastaEntry("File", "gi|155030192", 
				"MADNLSDTLKKLKITAVDKTEDSLEGCLDCLLQALAQNNTETSEKIQASGILQLFASLLTPQSSCKAKVA"+
				"NIIAEVAKNEFMRIPCVDAGLISPLVQLLNSKDQEVLLQTGRALGNICYDSHEGRSAVDQAGGAQIVIDH"+
				"LRSLCSITDPANEKLLTVFCGMLMNYSNENDSLQAQLINMGVIPTLVKLLGIHCQNAALTEMCLVAFGNL"+
				"AELESSKEQFASTNIAEELVKLFKKQIEHDKREMIFEVLAPLAENDAIKLQLVEAGLVECLLEIVQQKVD"+
				"SDKEDDITELKTGSDLMVLLLLGDESMQKLFEGGKGSVFQRVLSWIPSNNHQLQLAGALAIANFARNDAN"+
				"CIHMVDNGIVEKLMDLLDRHVEDGNVTVQHAALSALRNLAIPVINKAKMLSAGVTEAVLKFLKSEMPPVQ"+
				"FKLLGTLRMLIDAQAEAAEQLGKNVKLVERLVEWCEAKDHAGVMGESNRLLSALIRHSKSKDVIKTIVQS"+
				"GGIKHLVTMATSEHVIMQNEALVALALIAALELGTAEKDLESAKLVQILHRLLADERSAPEIKYNSMVLI"+
				"CALMGSECLHKEVQDLAFLDVVSKLRSHENKSVAQQASLTEQRLTVES");
		trie.addFasta(fasta);
		assertEquals(2, entries.get(0).getAccessions().size());
		assertEquals(1, entries.get(2).getAccessions().size());
	}
}
