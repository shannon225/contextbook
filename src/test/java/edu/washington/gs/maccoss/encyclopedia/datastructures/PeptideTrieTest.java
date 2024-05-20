package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Optional;
import java.util.StringTokenizer;

import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MSPReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PrositCSVWriter;
import edu.washington.gs.maccoss.encyclopedia.utils.io.LineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.LineParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import junit.framework.TestCase;

public class PeptideTrieTest extends TestCase {
	
	public static void main(String[] args) throws Exception {

		//File file=new File("/Users/searleb/Documents/encyclopedia/positional_isomers/phosphopedia_human_phosphopeptide_HCD_library.dlib");
		File file=new File("/Users/searleb/Downloads/human/pan_human_library.dlib");
		File fastaFile=new File("/Users/searleb/Downloads/uniprot_sprot.fasta");
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		LibraryFile library=new LibraryFile();
		library.openFile(file);
		
		final HashSet<String> motifs=new HashSet<String>();
		PeptideTrie<LibraryEntry> motifTrie=new PeptideTrie<LibraryEntry>(library.getAllEntries(false, params.getAAConstants())) {
			@Override
			protected void processMatch(FastaEntryInterface fasta, LibraryEntry entry, int start) {
				String peptideModSeq = entry.getPeptideModSeq();
				peptideModSeq=peptideModSeq.replace("T[+79.966331]", "T");
				peptideModSeq=peptideModSeq.replace("Y[+79.966331]", "Y");
				peptideModSeq=peptideModSeq.replace("[+15.994915]", "");
				peptideModSeq=peptideModSeq.replace("[+42.010565]", "");
				peptideModSeq=peptideModSeq.replace("[+57.0214635]", "");
				peptideModSeq=peptideModSeq.replace("[+57.021464]", "");
				peptideModSeq=peptideModSeq.replace("[+58.00548]", "");
				peptideModSeq=peptideModSeq.replace("[+99.0320285]", "");
				peptideModSeq=peptideModSeq.replace("[+121.976896]", "");
				
				if (peptideModSeq.indexOf('[')>=0) {
					System.out.println(peptideModSeq);
					System.exit(1);
				}
				
				int offset=peptideModSeq.indexOf("S");
				if (offset<0) return;
				
				int begin=start-5+offset;
				int end=start+5+offset+1;

				String sequence = fasta.getSequence();
				StringBuilder sb=new StringBuilder();
				for (int i = begin; i < end; i++) {
					if (i<0) {
						sb.append(' ');
					} else {
						if (i>=sequence.length()) {
							sb.append(' ');
						} else {
							sb.append(sequence.charAt(i));
						}
					}
				}
				motifs.add(sb.toString());
			}
		};

		ArrayList<FastaEntryInterface> proteins=FastaReader.readFasta(fastaFile, params);
		motifTrie.addFasta(proteins);
		
		for (String motif : motifs) {
			System.out.println(motif);
		}
	}

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
	
	public static void main2(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		float scoreThreshold = 1.00E-08f;
		File[] humanLibraries=new File[] {
				new File("/Users/searleb/Documents/damien/dda_library_search/HeLa.dlib"),
				new File("/Users/searleb/Documents/damien/dda_library_search/MCF7.dlib"),
				new File("/Users/searleb/Documents/damien/dda_library_search/pan_human_library.dlib"),
		};
		HashSet<PeptidePrecursor> expectedPeptides=new HashSet<>();
		for (File file : humanLibraries) {
			LibraryFile library=new LibraryFile();
			library.openFile(file);
			ArrayList<PeptidePrecursor> allPeptidePrecursors = library.getAllPeptidePrecursors(parameters.getAAConstants());
			for (PeptidePrecursor peptidePrecursor : allPeptidePrecursors) {
				if (peptidePrecursor.getPeptideSeq().length()<7||peptidePrecursor.getPeptideSeq().length()>30) continue;
				if (peptidePrecursor.getPrecursorCharge()>4) continue;
				
				expectedPeptides.add(new SimplePeptidePrecursor(peptidePrecursor.getPeptideSeq(), peptidePrecursor.getPrecursorCharge(), parameters.getAAConstants()));
			}
			System.out.println("New Expected Library: "+allPeptidePrecursors.size()+", Final Total "+expectedPeptides.size()+" added after "+file.getName());
		}
		System.out.println("Total Expected: "+expectedPeptides.size());
		
		File thegpmFile=new File("/Users/searleb/Documents/damien/dda_library_search/Homo_sapiens.tsv");
		File fastaFile = new File("/Users/searleb/Documents/damien/dda_library_search/uniprot_human-reference_reviewed_2022mar02.fasta");
		ArrayList<FastaEntryInterface> proteins=FastaReader.readFasta(fastaFile, parameters);
		
		int[] scores=new int[20];
		
		Collection<PeptidePrecursor> peptides=new HashSet<PeptidePrecursor>();
		LineParser.parseFile(thegpmFile, new LineParserMuscle() {
			boolean first=true;
			
			@Override
			public void processRow(String row) {
				if (first) {
					first=false;
					return;
				}
				StringTokenizer st=new StringTokenizer(row);
				String peptide=st.nextToken();
				if (peptide.length()<7||peptide.length()>30) return; // outside Prosit
				
				int p1Count=Integer.parseInt(st.nextToken());
				int p2Count=Integer.parseInt(st.nextToken());
				int p3Count=Integer.parseInt(st.nextToken());
				int p4Count=Integer.parseInt(st.nextToken());
				
				float score=Float.parseFloat(st.nextToken());
				if (score>scoreThreshold) return;
				
				int index=Math.round(-Log.protectedLog10(score));
				if (index<0) index=0;
				if (index>=scores.length) index=scores.length-1;
				scores[index]++;
				
				if (p1Count>0) peptides.add(new SimplePeptidePrecursor(peptide, (byte)1, parameters.getAAConstants()));
				if (p2Count>0) peptides.add(new SimplePeptidePrecursor(peptide, (byte)2, parameters.getAAConstants()));
				if (p3Count>0) peptides.add(new SimplePeptidePrecursor(peptide, (byte)3, parameters.getAAConstants()));
				if (p4Count>0) peptides.add(new SimplePeptidePrecursor(peptide, (byte)4, parameters.getAAConstants()));
			}
			
			@Override
			public void cleanup() {
			}
		});
		
		HashSet<String> uniqueSequences=new HashSet<String>();
		for (PeptidePrecursor peptide : peptides) {
			uniqueSequences.add(peptide.getPeptideSeq());
		}
		System.out.println("Starting: "+peptides.size()+" total precursors, "+uniqueSequences.size()+" unique sequences");
		
		int missingPrecursors=0;
		int missingPeptides=0;
		for (PeptidePrecursor pep : expectedPeptides) {
			if (!peptides.contains(pep)) {
				missingPrecursors++;
			}
			if (!uniqueSequences.contains(pep.getPeptideSeq())) {
				missingPeptides++;
			}
		}
		System.out.println("Missing precursors: "+missingPrecursors+", missing peptides: "+missingPeptides);
		
		System.out.println("Score distribution: ");
		for (int i = 0; i < scores.length; i++) {
			System.out.println(i+"\t"+scores[i]);
		}
		
		PeptidePrintingTrie trie=new PeptidePrintingTrie(peptides, Optional.ofNullable(DigestionEnzyme.getEnzyme("Trypsin")));
		trie.addFasta(proteins);
		
		uniqueSequences=new HashSet<String>();
		for (PeptidePrecursor peptide : trie.matchingPeptides) {
			uniqueSequences.add(peptide.getPeptideSeq());
		}
		System.out.println("Filtered: "+trie.matchingPeptides.size()+" total precursors, "+uniqueSequences.size()+" unique sequences");
		
		missingPrecursors=0;
		missingPeptides=0;
		for (PeptidePrecursor pep : expectedPeptides) {
			if (!trie.matchingPeptides.contains(pep)) {
				missingPrecursors++;
			}
			if (!uniqueSequences.contains(pep.getPeptideSeq())) {
				missingPeptides++;
			}
		}
		System.out.println("Filtered Missing precursors: "+missingPrecursors+", missing peptides: "+missingPeptides);
		//System.out.println(scoreThreshold+"\t"+trie.matchingPeptides.size()+"\t"+missingPrecursors/(float)expectedPeptides.size());
		
		HashSet<PeptidePrecursor> precursors=new HashSet<>(trie.matchingPeptides);
		precursors.addAll(expectedPeptides);
		System.out.println("Final library size: "+precursors.size());
		for (int nce : new int[] {17, 26, 29, 33, 35, 38, 46, 50}) {
			String prositOutputFileName = "/Users/searleb/Documents/damien/dda_library_search/prosit_input_nce_"+nce+".csv";
			PrositCSVWriter.writePrositFile(prositOutputFileName, nce, (byte)2, false, false, precursors);
		}
		
	}
	
	private static class PeptidePrintingTrie extends PeptideTrie<PeptidePrecursor> {
		private final HashSet<PeptidePrecursor> matchingPeptides=new HashSet<PeptidePrecursor>();
		
		
		public PeptidePrintingTrie(Collection<PeptidePrecursor> entries, Optional<DigestionEnzyme> enzyme) {
			super(entries, enzyme);
		}


		@Override
		protected void processMatch(FastaEntryInterface fasta, PeptidePrecursor entry, int start) {
			matchingPeptides.add(entry);
		}
	}
}
