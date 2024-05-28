package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SimplePeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PrositCSVWriter;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import gnu.trove.map.hash.TObjectIntHashMap;

public class TheGPMParserTest {

	public static void main(String[] args) throws Exception {
		//File gpmFile=new File("/Users/searleb/Documents/encyclopedia/thegpm_peptide_lists/Eukaryotes/Homo_sapiens.tsv");
		//File fasta=new File("/Users/searleb/Documents/encyclopedia/thegpm_peptide_lists/fastas/human.fasta");
		SearchParameters params = SearchParameterParser.getDefaultParametersObject();
		
		List<DigestionEnzyme> enzymes=new ArrayList<DigestionEnzyme>();
		DigestionEnzyme trypsin = DigestionEnzyme.getEnzyme("Trypsin");
		enzymes.add(trypsin);
		enzymes.add(DigestionEnzyme.getEnzyme("Lys-C"));
		enzymes.add(DigestionEnzyme.getEnzyme("Lys-N"));
		enzymes.add(DigestionEnzyme.getEnzyme("Arg-C"));
		enzymes.add(DigestionEnzyme.getEnzyme("Glu-C"));
		enzymes.add(DigestionEnzyme.getEnzyme("Chymotrypsin"));
		enzymes.add(DigestionEnzyme.getEnzyme("Pepsin A"));
		enzymes.add(DigestionEnzyme.getEnzyme("Elastase"));
		enzymes.add(DigestionEnzyme.getEnzyme("Thermolysin"));
		
		HashMap<DigestionEnzyme, HashSet<GPMPeptide>> peptideSetByEnzyme=new HashMap<>();
		for (DigestionEnzyme enzyme : enzymes) {
			peptideSetByEnzyme.put(enzyme, new HashSet<>());
		}

		File dir=new File("/Users/searleb/Documents/encyclopedia/thegpm_peptide_lists/");
		File gpmFile=new File(dir, "cRAP.tsv");
		File fasta=new File(dir, "crap.fasta");
		
		HashMap<String, GPMPeptide> pepBySeq = getPeptideSet(params, gpmFile);
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(fasta, params);
		
		for (DigestionEnzyme enzyme : enzymes) {
			HashSet<GPMPeptide> kept = matchToEnzyme(entries, params, pepBySeq, enzyme, 999999);
			peptideSetByEnzyme.get(enzyme).addAll(kept);
			System.out.println(enzyme+" --> "+peptideSetByEnzyme.get(enzyme).size());
		}
		
		gpmFile=new File(dir, "Eukaryotes/Homo_sapiens.tsv");
		fasta=new File(dir, "fastas/human.fasta");
		
		pepBySeq = getPeptideSet(params, gpmFile);
		entries=FastaReader.readFasta(fasta, params);
		
		for (DigestionEnzyme enzyme : enzymes) {
			HashSet<GPMPeptide> kept = matchToEnzyme(entries, params, pepBySeq, enzyme, 2500);
			peptideSetByEnzyme.get(enzyme).addAll(kept);
			System.out.println(enzyme+" --> "+peptideSetByEnzyme.get(enzyme).size());
		}
		
		gpmFile=new File(dir, "Eukaryotes/Mus_musculus.tsv");
		fasta=new File(dir, "fastas/mouse.fasta");
		
		pepBySeq = getPeptideSet(params, gpmFile);
		entries=FastaReader.readFasta(fasta, params);
		
		for (DigestionEnzyme enzyme : enzymes) {
			HashSet<GPMPeptide> kept = matchToEnzyme(entries, params, pepBySeq, enzyme, 2500);
			peptideSetByEnzyme.get(enzyme).addAll(kept);
			System.out.println(enzyme+" --> "+peptideSetByEnzyme.get(enzyme).size());
		}
		
		gpmFile=new File(dir, "Eukaryotes/Saccharomyces_cerevisiae.tsv");
		fasta=new File(dir, "fastas/yeast.fasta");
		
		pepBySeq = getPeptideSet(params, gpmFile);
		entries=FastaReader.readFasta(fasta, params);
		
		for (DigestionEnzyme enzyme : enzymes) {
			HashSet<GPMPeptide> kept = matchToEnzyme(entries, params, pepBySeq, enzyme, 2500);
			peptideSetByEnzyme.get(enzyme).addAll(kept);
			System.out.println(enzyme+" --> "+peptideSetByEnzyme.get(enzyme).size());
		}
		
		gpmFile=new File(dir, "Eukaryotes/Rattus_norvegicus.tsv");
		fasta=new File(dir, "fastas/rat.fasta");
		
		pepBySeq = getPeptideSet(params, gpmFile);
		entries=FastaReader.readFasta(fasta, params);
		
		for (DigestionEnzyme enzyme : enzymes) {
			HashSet<GPMPeptide> kept = matchToEnzyme(entries, params, pepBySeq, enzyme, 1000);
			peptideSetByEnzyme.get(enzyme).addAll(kept);
			System.out.println(enzyme+" --> "+peptideSetByEnzyme.get(enzyme).size());
		}
		
		gpmFile=new File(dir, "Eukaryotes/Drosophila_melanogaster.tsv");
		fasta=new File(dir, "fastas/fruitfly.fasta");
		
		pepBySeq = getPeptideSet(params, gpmFile);
		entries=FastaReader.readFasta(fasta, params);
		
		for (DigestionEnzyme enzyme : enzymes) {
			HashSet<GPMPeptide> kept = matchToEnzyme(entries, params, pepBySeq, enzyme, 1000);
			peptideSetByEnzyme.get(enzyme).addAll(kept);
			System.out.println(enzyme+" --> "+peptideSetByEnzyme.get(enzyme).size());
		}
		
		gpmFile=new File(dir, "Prokaryotes/Escherichia_coli_K-12.tsv");
		fasta=new File(dir, "fastas/ecoli.fasta");
		
		pepBySeq = getPeptideSet(params, gpmFile);
		entries=FastaReader.readFasta(fasta, params);
		
		for (DigestionEnzyme enzyme : enzymes) {
			HashSet<GPMPeptide> kept = matchToEnzyme(entries, params, pepBySeq, enzyme, 1000);
			peptideSetByEnzyme.get(enzyme).addAll(kept);
			System.out.println(enzyme+" --> "+peptideSetByEnzyme.get(enzyme).size());
		}
		

		for (DigestionEnzyme enzyme : enzymes) {
			HashSet<PeptidePrecursor> kept=new HashSet<>(peptideSetByEnzyme.get(enzyme));
			
			File f=new File(dir, "prosit/"+enzyme.getName()+".csv");
			PrositCSVWriter.writePrositFile(f.getAbsolutePath(), 33, (byte)3, false, false, kept);
		}
	}

	private static HashMap<String, GPMPeptide> getPeptideSet(SearchParameters params, File gpmFile) {
		HashSet<GPMPeptide> set=getPeptidesFromGPMList(gpmFile, params.getAAConstants());
		HashMap<String, GPMPeptide> pepBySeq=new HashMap<>();
		for (GPMPeptide pep : set) {
			pepBySeq.put(pep.getPeptideSeq(), pep);
		}
		return pepBySeq;
	}

	private static HashSet<GPMPeptide> matchToEnzyme(ArrayList<FastaEntryInterface> entries, SearchParameters params,
			HashMap<String, GPMPeptide> pepBySeq, DigestionEnzyme enzyme, int topPeptideKept) {
		HashSet<GPMPeptide> matching=new HashSet<>();
		TObjectIntHashMap<FastaEntryInterface> proteinCounter=new TObjectIntHashMap<>();
		for (FastaEntryInterface entry : entries) {
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 7, 30, 1, params.getAAConstants(), false);
			for (FastaPeptideEntry peptide : peptides) {
				GPMPeptide pep=pepBySeq.get(peptide.getSequence());
				if (pep!=null) {
					matching.add(pep);
					int count=pep.getFrequency();
					proteinCounter.adjustOrPutValue(entry, count, count);
				}
			}
		}

		ArrayList<GPMPeptide> peps=new ArrayList<>(matching);
		Collections.sort(peps);
		Collections.reverse(peps);

		HashSet<GPMPeptide> kept=new HashSet<>();
		for (int i = 0; i < Math.min(topPeptideKept, peps.size()); i++) {
			kept.add(peps.get(i));
		}

//		ArrayList<ScoredObject<String>> scores=new ArrayList<>();
//		proteinCounter.forEachEntry(new TObjectIntProcedure<FastaEntryInterface>() {
//			@Override
//			public boolean execute(FastaEntryInterface a, int b) {
//				scores.add(new ScoredObject<String>(b, a.getAccession()));
//				return true;
//			}
//		});
//		
//		Collections.sort(scores);
//		Collections.reverse(scores);
		
//		System.out.println(enzyme.toString()+": "+peps.size()+" total peptides");
//		for (int i = 0; i < 10; i++) {
//			System.out.println("\t"+scores.get(i).toString());
//		}
//
//		for (int i = 0; i < 10; i++) {
//			GPMPeptide gpmPeptide = peps.get(i);
//			System.out.println("\t"+gpmPeptide.getFrequency()+", "+gpmPeptide.toString());
//		}
//		System.out.println();
		return kept;
	}

	private static HashSet<GPMPeptide> getPeptidesFromGPMList(File f, AminoAcidConstants aaConstants) {
		final HashSet<GPMPeptide> set=new HashSet<TheGPMParserTest.GPMPeptide>();
		
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String sequence=row.get("sequence");
				int z1=Integer.parseInt(row.get("z=1"));
				int z2=Integer.parseInt(row.get("z=2"));
				int z3=Integer.parseInt(row.get("z=3"));
				int z4=Integer.parseInt(row.get("z=4"));
				float e=Float.parseFloat(row.get("E"));
				set.add(new GPMPeptide(sequence, z1, z2, z3, z4, e, aaConstants));
			}
			
			@Override
			public void cleanup() {
			}
		};
		
		TableParser.parseTSV(f, muscle);
		
		return set;
	}
	
	private static class GPMPeptide extends SimplePeptidePrecursor {
		private final int z1;
		private final int z2;
		private final int z3;
		private final int z4;
		private final float e;
		
		public GPMPeptide(String sequence, int z1, int z2, int z3, int z4, float e, AminoAcidConstants aaConstants) {
			super(sequence, getBestPrecursorCharge(z1, z3, z4, z4), aaConstants);
			this.z1 = z1;
			this.z2 = z2;
			this.z3 = z3;
			this.z4 = z4;
			this.e = e;
		}
		
		public int getFrequency() {
			switch (getPrecursorCharge()) {
				case 1: return z1;
				case 2: return z2;
				case 3: return z3;
				case 4: return z4;
			}
			return 0;
		}
		
		@Override
		public int compareTo(PeptidePrecursor o) {
			if (o instanceof GPMPeptide) {
				int c=Integer.compare(getFrequency(), ((GPMPeptide)o).getFrequency());
				if (c!=0) return c;
			}
			return super.compareTo(o);
		}
		
		@Override
		public String toString() {
			return getPeptideSeq()+": "+z1+","+z2+","+z3+","+z4+" --> "+e;
		}
		
		public static byte getBestPrecursorCharge(int z1, int z2, int z3, int z4) {
			int max=z1;
			byte charge=1;
			if (max<z2) {
				max=z2;
				charge=2;
			}
			if (max<z3) {
				max=z3;
				charge=3;
			}
			if (max<z4) {
				max=z4;
				charge=4;
			}
			return charge;
		}
	}
}
