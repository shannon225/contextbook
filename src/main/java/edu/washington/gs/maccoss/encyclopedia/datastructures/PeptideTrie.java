package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Collection;

import gnu.trove.map.hash.TCharObjectHashMap;

public class PeptideTrie {
	TrieNode head=new TrieNode('$');
	
	/**
	 * peptide trie stored backwards (so K/R comes first)
	 * @param entries
	 */
	public PeptideTrie(Collection<LibraryEntry> entries) {
		for (LibraryEntry entry : entries) {
			char[] sequence=entry.getPeptideSeq().toCharArray();
			
			TrieNode node=head;
			for (int i=sequence.length-1; i>=0; i--) {
				node=node.getOrCreate(sequence[i]);
			}
			node.addEntry(entry);
		}
	}
	
	public void addFasta(FastaEntryInterface fasta) {
		String accession=fasta.getAccession();
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
					for (LibraryEntry entry : node.entries) {
						entry.getAccessions().add(accession);
					}
				}
			}
		}
	}

	private class TrieNode {
		private final char aa;
		private final TCharObjectHashMap<TrieNode> children=new TCharObjectHashMap<PeptideTrie.TrieNode>();
		private final ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();

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
		
		public void addEntry(LibraryEntry entry) {
			entries.add(entry);
		}
	}
}
