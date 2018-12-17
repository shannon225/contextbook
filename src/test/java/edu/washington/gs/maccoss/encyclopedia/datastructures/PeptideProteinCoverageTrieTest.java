package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.map.hash.TIntFloatHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;

public class PeptideProteinCoverageTrieTest {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50),
			DigestionEnzyme.getEnzyme("trypsin"), false, true, false);
	
	public static void main(String[] args) throws Exception {
		File fastaFile=new File("/Users/searleb/Downloads/IARPA_20181214_10_pepFDR.fasta");
		File chromlib=new File("/Users/searleb/Downloads/10pFDR.elib");
		LibraryFile lib=new LibraryFile();
		lib.openFile(chromlib);
		
		TObjectIntHashMap<String> groupByAccession=new TObjectIntHashMap<>();
		try (Connection c=lib.getConnection()) {
			try (PreparedStatement s=c.prepareStatement("select proteingroup, proteinaccession from proteinscores where isdecoy=0")) {

				ResultSet rs=s.executeQuery();
				while (rs.next()) {
					int group=rs.getInt(1);
					String accession=rs.getString(2);
					groupByAccession.put(accession, group);
				}
			}
		}
		System.out.println("Found "+groupByAccession.size()+" accessions");

		ArrayList<LibraryEntry> allEntries=lib.getAllEntries(false, PARAMETERS.getAAConstants());
		ArrayList<FastaEntryInterface> fasta=FastaReader.readFasta(fastaFile, PARAMETERS);
		PeptideAccessionMatchingTrie fastaLookup=new PeptideAccessionMatchingTrie(allEntries);
		fastaLookup.addFasta(fasta);
		
		PeptideProteinCoverageTrie trie=new PeptideProteinCoverageTrie(allEntries);
		trie.addFasta(fasta);

		TIntFloatHashMap bestCoverageForGroup=new TIntFloatHashMap();
		TIntObjectHashMap<String> bestAccessionForGroup=new TIntObjectHashMap<>();
		HashMap<String, boolean[]> map=trie.getCoverageMap();
		for (Entry<String, boolean[]> entry : map.entrySet()) {
			String accession=entry.getKey();
			boolean[] value=entry.getValue();
			int trues=0;
			for (int i=0; i<value.length; i++) {
				if (value[i]) trues++;
			}
			float coverage=trues/(float)value.length;
			int group=groupByAccession.get(accession);
			if (group>0&&coverage>bestCoverageForGroup.get(group)) { // will return 0 if missing
				bestCoverageForGroup.put(group, coverage);
				bestAccessionForGroup.put(group, accession);
			}
		}
		System.out.println("Found "+bestAccessionForGroup.size()+" groups");
		
		for (String accession : bestAccessionForGroup.valueCollection()) {
			boolean[] value=map.get(accession);
			int trues=0;
			for (int i=0; i<value.length; i++) {
				if (value[i]) trues++;
			}
			float coverage=trues/(float)value.length;

			System.out.println(accession+"\t"+trie.getPeptideCount().get(accession)+"\t"+trie.getUniquePeptideCount().get(accession)+"\t"+value.length+"\t"+coverage);
		}
	}
}
