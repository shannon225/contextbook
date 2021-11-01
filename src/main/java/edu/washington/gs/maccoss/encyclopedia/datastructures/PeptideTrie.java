package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import gnu.trove.map.hash.TCharObjectHashMap;

public abstract class PeptideTrie <T extends PeptidePrecursor> {
	TrieNode head=new TrieNode('$');
	private final Optional<DigestionEnzyme> enzyme;
	
	public PeptideTrie(Collection<T> entries) {
		this(entries, Optional.empty());
	}
	
	/**
	 * peptide trie stored backwards (so K/R comes first)
	 * @param entries
	 */
	public PeptideTrie(Collection<T> entries, Optional<DigestionEnzyme> enzyme) {
		this.enzyme=enzyme;
		for (T entry : entries) {
			char[] sequence=entry.getPeptideSeq().toCharArray();
			
			TrieNode node=head;
			for (int i=sequence.length-1; i>=0; i--) {
				node=node.getOrCreate(sequence[i]);
			}
			node.addEntry(entry);
		}
	}
	
	public void addFasta(ArrayList<FastaEntryInterface> fasta) {
		for (FastaEntryInterface f : fasta) {
			addFasta(f);
		}
	}
	
	public void addFasta(FastaEntryInterface fasta) {
		char[] sequence=fasta.getSequence().toCharArray();
		for (int i=sequence.length-1; i>=0; i--) {
			// move trie along sequence
			TrieNode node=head.get(sequence[i]);
			if (node==null) continue;

			// process current reverse peptide
			for (int j=i-1; j>=0; j--) {
				node=node.get(sequence[j]);
				if (node==null) break;
				
				if (node.entries.size()>0) {
					if (j==0||!enzyme.isPresent()||enzyme.get().isCutSite(sequence[j-1], sequence[j])) {
						// either the beginning of the protein, the enzyme isn't used, or the enzyme indicates a cut site
						for (T entry : node.entries) {
							// add fasta protein entry to each peptide entry
							processMatch(fasta, entry, j);
						}
					} else if (enzyme.isPresent()) {
						System.out.println(sequence[j-1]+","+sequence[j]+"="+node.entries.get(0).getPeptideSeq().charAt(0)+" "+fasta.getAccession());
					}
				}
			}
		}
	}

	protected abstract void processMatch(FastaEntryInterface fasta, T entry, int start);

	private class TrieNode {
		private final char aa;
		private final TCharObjectHashMap<TrieNode> children=new TCharObjectHashMap<TrieNode>();
		private final ArrayList<T> entries=new ArrayList<T>();

		public TrieNode(char aa) {
			this.aa=aa;
		}
		@Override
		public String toString() {
			return Character.toString(aa);
		}
		
		/**
		 * not null protected!
		 */
		public TrieNode get(char aa) {
			return children.get(aa);
		}
		
		public TrieNode getOrCreate(char aa) {
			TrieNode node=get(aa);
			if (node==null) {
				node=new TrieNode(aa);
				children.put(aa, node);
			}
			return node;
		}
		
		public void addEntry(T entry) {
			entries.add(entry);
		}
	}
}
