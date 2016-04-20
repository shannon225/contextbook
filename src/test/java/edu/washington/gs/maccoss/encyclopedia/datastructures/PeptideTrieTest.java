package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.InputStream;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.filereaders.MSPReader;
import junit.framework.TestCase;

public class PeptideTrieTest extends TestCase {

	public void testTrie() {
		InputStream is=getClass().getResourceAsStream("/truncated.msp");
		ArrayList<LibraryEntry> entries=MSPReader.readMSP(is, "truncated.msp");
		PeptideTrie trie=new PeptideTrie(entries);
		
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
