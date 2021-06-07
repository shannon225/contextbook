package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.zip.DataFormatException;

import org.relaxng.datatype.DatatypeException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideAccessionMatchingTrie;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class LibraryEntryCleaner {
	public static void main(String[] args) throws Exception {
//		File originalLibraryFile=new File("/Users/searleb/Documents/iarpa/new_llnl_individual/IARPA_3clib_plus_llnl.trypsin.z3_nce33.dlib");
//		File newLibraryFile=new File("/Users/searleb/Documents/iarpa/new_llnl_individual/uniq_IARPA_3clib_plus_llnl.trypsin.z3_nce33.dlib");
//		File fastaFile=new File("/Users/searleb/Documents/iarpa/new_llnl_individual/IARPA_20201215_plus_ref_vars.fasta");
//		
//		cleanLibrary(true, originalLibraryFile, newLibraryFile, fastaFile, SearchParameterParser.getDefaultParametersObject());
		

		//File fileName=new File("/Users/searleb/Documents/grants/phosphotau mutagenesis/barcodes_prosit.dlib");
		//File newFileName=new File("/Users/searleb/Documents/grants/phosphotau mutagenesis/ecoli_barcodes_prosit.dlib");
		//File fasta=new File("/Users/searleb/Documents/grants/phosphotau mutagenesis/uniprot_ecoli_27dec2020.fasta");

		//fileName = new File("/Users/searleb/Documents/cobbs/2020dec30_cobbs/2020dec03_cobbs_cmv_inf_gpfdia_clib.elib");
		//newFileName = new File("/Users/searleb/Documents/cobbs/2020dec30_cobbs/2020dec03_cobbs_cmv_inf_gpfdia_clib_hcmv_only.elib");
		//fasta=new File("/Users/searleb/Documents/cobbs/2020dec30_cobbs/merlin.fasta");
		
		File fileName=new File("/Users/searleb/Documents/iarpa/new_llnl_individual_2/unclean_IARPA_20210401_3clib_plus_llnl.trypsin.z3_nce33.dlib");
		File newFileName=new File("/Users/searleb/Documents/iarpa/new_llnl_individual_2/clean_IARPA_20210401_3clib_plus_llnl.trypsin.z3_nce33.dlib");
		File fasta=new File("/Users/searleb/Documents/iarpa/new_llnl_individual_2/IARPA_20210401.fasta_plus_llnl_ref_vars.fasta");
		
		cleanLibrary(true, fileName, newFileName, fasta, SearchParameterParser.getDefaultParametersObject());
	}
	
	public static boolean doesLibraryNeedCleaning(LibraryFile libraryFile) throws IOException, SQLException {
		Connection c = libraryFile.getConnection();
		try {
			Statement s=c.createStatement();
			try {
				Logger.logLine("Getting count of protein accessions...");
				ResultSet rs=s.executeQuery("select count(ProteinAccession) from peptidetoprotein");
				int count=rs.getInt(0);
				Logger.logLine("...found "+count);
				
				if (count>1) {
					Logger.logLine(libraryFile.getName()+" does not need to be cleaned!");
					return false;
				}
				rs.close();

			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
		Logger.logLine(libraryFile.getName()+" needs to be cleaned!");
		return true;
	}
	
	public static LibraryFile cleanLibrary(boolean smallerScoresAreGood, File originalLibraryFile, File newLibraryFile, File fastaFile, SearchParameters parameters) throws IOException, DatatypeException, DataFormatException, SQLException {
		LibraryFile oldLibrary=new LibraryFile();
		oldLibrary.openFile(originalLibraryFile);
		ArrayList<LibraryEntry> originalEntries=oldLibrary.getAllEntries(false, parameters.getAAConstants());
		
		HashMap<String, LibraryEntry> uniqueEntryMap=new HashMap<>();
		for (LibraryEntry entry : originalEntries) {
			entry.getAccessions().clear(); // clear accessions so they can be re-populated
			String key=entry.getPeptideModSeq()+"_"+entry.getPrecursorCharge();
			
			LibraryEntry prev=uniqueEntryMap.get(key);
			if (prev==null||(smallerScoresAreGood?(prev.getScore()>entry.getScore()):(prev.getScore()<entry.getScore()))) {
				uniqueEntryMap.put(key, entry);
			}
		}
		ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>(uniqueEntryMap.values());
		
		Logger.logLine("Found "+entries.size()+" unique peptide entries");

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
		
		LibraryFile newLibrary=new LibraryFile();
		newLibrary.openFile();
		Logger.logLine("Writing library file "+newLibrary.getName());
		newLibrary.dropIndices();
		newLibrary.addEntries(entries);
		newLibrary.addProteinsFromEntries(entries);
		newLibrary.createIndices();
		newLibrary.saveAsFile(newLibraryFile);
		return newLibrary;
	}
}
