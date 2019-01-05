package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.Collection;
import java.util.HashMap;

import gnu.trove.map.hash.TObjectIntHashMap;

public class PeptideProteinCoverageTrie extends PeptideTrie<LibraryEntry> {
	HashMap<String, boolean[]> coverageMap=new HashMap<>();
	TObjectIntHashMap<String> peptideCount=new TObjectIntHashMap<>();
	TObjectIntHashMap<String> uniquePeptideCount=new TObjectIntHashMap<>();
	
	public PeptideProteinCoverageTrie(Collection<LibraryEntry> entries) {
		super(entries);
	}

	@Override
	protected void processMatch(FastaEntryInterface fasta, LibraryEntry entry, int start) {
		peptideCount.adjustOrPutValue(fasta.getAccession(), 1, 1);
		if (entry.getAccessions().size()==1) {
			uniquePeptideCount.adjustOrPutValue(fasta.getAccession(), 1, 1);
		}
		boolean[] map=coverageMap.get(fasta.getAccession());
		if (map==null) {
			map=new boolean[fasta.getSequence().length()];
			coverageMap.put(fasta.getAccession(), map);
		}
		for (int i=start; i<start+entry.getPeptideSeq().length(); i++) {
			if (i>=0&&i<map.length) {
				map[i]=true;
			}
		}
	}
	
	public HashMap<String, boolean[]> getCoverageMap() {
		return coverageMap;
	}
	
	public TObjectIntHashMap<String> getPeptideCount() {
		return peptideCount;
	}
	
	public TObjectIntHashMap<String> getUniquePeptideCount() {
		return uniquePeptideCount;
	}
}
