package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filewriters.FastaWriter;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import gnu.trove.list.array.TCharArrayList;

public class FastaDatabaseReverser {

	public static void main(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		File f=new File("/Users/searleb/Documents/school/uniprot-9606.fasta");
		FastaWriter revwriter=new FastaWriter(new File(f.getAbsolutePath()+".rev.fasta"));
		FastaWriter revprowriter=new FastaWriter(new File(f.getAbsolutePath()+".revpro.fasta"));
		FastaWriter randprowriter=new FastaWriter(new File(f.getAbsolutePath()+".randpro.fasta"));
		FastaWriter randdbwriter=new FastaWriter(new File(f.getAbsolutePath()+".randdb.fasta"));
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);

		int[] aaCounts=new int[26];
		for (FastaEntryInterface entry : entries) {
			String sequence=entry.getSequence();
			for (char c : sequence.toCharArray()) {
				aaCounts[getIndex(c)]++;
			}
		}
		int aaTotalCounts=General.sum(aaCounts);
		
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		int seed=7*7*7*7*7;
		for (FastaEntryInterface entry : entries) {
			seed=RandomGenerator.randomInt(seed);
			char[] sequence=entry.getSequence().toCharArray();
			StringBuilder newSequence=new StringBuilder();
			
			TCharArrayList buffer=new TCharArrayList();
			for (int i = 0; i < sequence.length; i++) {
				if (enzyme.isTargetPreSite(sequence[i])) {
					newSequence.append(shuffle(buffer.toArray(), seed));
					seed=RandomGenerator.randomInt(seed);
					newSequence.append(sequence[i]);
					buffer.clear();
				} else {
					buffer.add(sequence[i]);
				}
			}
			if (buffer.size()>0) {
				buffer.reverse();
				newSequence.append(buffer.toArray());
			}

			revwriter.write(new FastaEntry(null, "RAND_"+entry.getAccession(), newSequence.toString()));
			revprowriter.write(new FastaEntry(null, "REVPRO_"+entry.getAccession(), new String(reverse(entry.getSequence().toCharArray()))));
			randprowriter.write(new FastaEntry(null, "RANDPRO_"+entry.getAccession(), new String(shuffle(entry.getSequence().toCharArray(), seed))));
			
			char[] random=new char[sequence.length];
			for (int i = 0; i < random.length; i++) {
				seed=RandomGenerator.randomInt(seed);
				random[i]=getRandomAA(aaCounts, aaTotalCounts, RandomGenerator.random(seed));
			}
			randdbwriter.write(new FastaEntry(null, "RANDDB_"+entry.getAccession(), new String(random)));
		}
		
		revwriter.close();
		revprowriter.close();
		randprowriter.close();
		randdbwriter.close();
	}
	
	private static char getRandomAA(int[] aaCounts, int aaTotalCounts, float prob) {
		int sum=(int)(prob*aaTotalCounts);
		
		int runningTotal=0;
		for (int i = 0; i < aaCounts.length; i++) {
			runningTotal+=aaCounts[i];
			if (runningTotal>sum) {
				return getChar(i);
			}
		}
		return getChar(aaCounts.length-1);
	}

	/**
	 * damages array!
	 * @param aas
	 * @return
	 */
	private static char[] reverse(char[] aas) {
		for (int i = 0, j = aas.length - 1; i < j; i++, j--) {
			char tmp = aas[i];
			aas[i] = aas[j];
			aas[j] = tmp;
		}
		return aas;
	}

	/**
	 * damages array!
	 * @param aas
	 * @return
	 */
	public static char[] shuffle(char[] aas, int shuffleSeed) {
		if (aas.length<=1) return aas;
		int start=0;
		int stop=aas.length-1;
		int diff=(stop)-start;
		
		// String.hashCode() is cross-platform consistent
		int seed=RandomGenerator.randomInt(shuffleSeed)+new String(aas).hashCode();
		for (int i=0; i<aas.length; i++) {
			seed=RandomGenerator.randomInt(seed);
			int index1=start+Math.abs(seed%diff+1);
			
			seed=RandomGenerator.randomInt(seed);
			int index2=start+Math.abs(seed%diff+1);
			if (index1!=index2) {
				char c=aas[index1];
				aas[index1]=aas[index2];
				aas[index2]=c;
			}
		}
		return aas;
	}
	
	private static int getIndex(char c) {
		assert(c>=65&&c<=90);
		return c-65;
	}
	
	private static char getChar(int index) {
		return (char)(index+65);
	}
}
