package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideAccessionMatchingTrie;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class TraMLSAXToLibraryConverter {
	public static void main(String[] args) {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		File tramlFile=new File("/Users/searleb/Documents/encyclopedia/bugs/traml/psgs_standard_consensus_filtered.TraML");
		File fastaFile=new File("/Users/searleb/Documents/encyclopedia/bugs/traml/PSGS_reference_peptides.fasta");
		Long time=System.currentTimeMillis();
		convertTraML(tramlFile, fastaFile, params);
		System.out.println(System.currentTimeMillis()-time);
	}

	public static LibraryFile convertTraML(File tramlFile, File fastaFile, SearchParameters parameters) {
		String absolutePath=tramlFile.getAbsolutePath();
		File libraryFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		return convertTraML(tramlFile, fastaFile, libraryFile, parameters);
	}
	
	public static LibraryFile convertTraML(File tramlFile, File fastaFile, File libraryFile, SearchParameters parameters) {
		try {
			TraMLSAXToLibraryProducer reader=new TraMLSAXToLibraryProducer(tramlFile, parameters);
			reader.run();
			ArrayList<LibraryEntry> entries=reader.getEntries();
			
			if (reader.hadError()) {
				Logger.errorLine("Unexpected parsing error when reading TraML!");
				throw new EncyclopediaException(reader.getError());
			}

			if (fastaFile!=null) {
				Logger.logLine("Reading Fasta file "+fastaFile.getName());
				ArrayList<FastaEntryInterface> proteins=FastaReader.readFasta(fastaFile, parameters);
			
				Logger.logLine("Constructing trie from library peptides");
				PeptideAccessionMatchingTrie trie=new PeptideAccessionMatchingTrie(entries);
				trie.addFasta(proteins);
			}

			int[] counts=new int[21];
			for (LibraryEntry entry : entries) {
				int size=Math.min(counts.length-1, entry.getAccessions().size());
				counts[size]++;
			}
			Logger.logLine("Accession count histogram: ");
			for (int i=0; i<counts.length; i++) {
				Logger.logLine(i+" Acc\t"+counts[i]+" Counts");
			}

			if (counts[0]>0) {
				Logger.errorLine(counts[0]+" library entries can't be linked to proteins! These entries will be dropped.");
			}
			
			LibraryFile library = new LibraryFile();
			library.openFile();
			Logger.logLine("Writing library file "+library.getName());
			library.dropIndices();
			library.addEntries(entries);
			library.addProteinsFromEntries(entries);
			library.createIndices();
			library.saveAsFile(libraryFile);
			return library;
			
		} catch (IOException e) {
			Logger.errorLine("IO Error parsing TraML:");
			Logger.errorException(e);
			throw new EncyclopediaException(e);
		} catch (SQLException e) {
			Logger.errorLine("SQL Error parsing TraML:");
			Logger.errorException(e);
			throw new EncyclopediaException(e);
		}
	}
}
