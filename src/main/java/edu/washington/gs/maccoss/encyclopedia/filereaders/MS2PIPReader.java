package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideAccessionMatchingTrie;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SimplePeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.filewriters.MS2PIPWriter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;

/**
 * Converter class for producing DLIB files from an MS2PIP
 * PEPREC and CSV file. See {@link MS2PIPWriter} for a class
 * that can generate PEPREC files from FASTA or library files.
 */
public class MS2PIPReader {
	private static final AminoAcidConstants emptyAAConstants=new AminoAcidConstants();

	public static LibraryFile convertMS2PIP(File peprecFile, File csvReportFile, File fastaFile, SearchParameters parameters) {
		String absolutePath=csvReportFile.getAbsolutePath();
		File libraryFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		
		HashMap<String, SimplePeptidePrecursor> precursorCodeMap=readInputFile(peprecFile);
		

		try {
			return convertFromMS2PIP(precursorCodeMap, csvReportFile, fastaFile, libraryFile, parameters);
		} catch (Exception e) {
			Logger.errorLine("Unexpected error parsing MS2PIP CSV:");
			Logger.errorException(e);
			throw new EncyclopediaException(e);
		}
	}
	
	public static LibraryFile convertFromMS2PIP(HashMap<String, SimplePeptidePrecursor> precursorCodeMap, File csvReportFile, File fastaFile, File libraryFile, SearchParameters parameters) throws IOException, SQLException {
		HashMap<String, PeptideEntry> entriesByCode=new HashMap<>();

		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String id=row.get("spec_id");
				PeptideEntry entry=entriesByCode.get(id);
				if (entry==null) {
					SimplePeptidePrecursor precursor=precursorCodeMap.get(id);
					if (precursor==null) {
						Logger.errorLine("Found missing peptide id: ["+id+"]");
						return;
					}
					float rt=Float.parseFloat(row.get("rt"));
					entry=new PeptideEntry(precursor.getPeptideModSeq(), precursor.getPrecursorCharge(), rt, csvReportFile.getName());
					entriesByCode.put(id, entry);
				}
				
				double mz=Double.parseDouble(row.get("mz")); 
				if (mz==0) {
					// MS2PIP BUG: MS2PIP produces +2H ions for +1H precursors with non-zero intensities. These obviously aren't in our fragmentation model. 
					//Logger.errorLine("Found unexpected ion type: ["+indexedIon+"] for peptide ["+id+", "+entry.peptideModSeq+"]");
				} else {
					float intensity=Float.parseFloat(row.get("prediction"));
					entry.addPeak(new Peak(mz, intensity));
				}
			}

			@Override
			public void cleanup() {
			}
		};
		TableParser.parseCSV(csvReportFile, muscle);

		ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
		for (PeptideEntry pe : entriesByCode.values()) {
			double precursorMZ=emptyAAConstants.getChargedMass(pe.peptideModSeq, pe.charge);
			HashSet<String> accessions=new HashSet<>();
			
			if (fastaFile==null) {
				accessions.add(PeptideUtils.getPeptideSeq(pe.peptideModSeq));
			}
			
			ImmutablePeptideEntry peptide=new ImmutablePeptideEntry(pe);
			LibraryEntry entry=new LibraryEntry(peptide.sourceFile, accessions, precursorMZ, peptide.charge, peptide.peptideModSeq, 1, peptide.rt, 0.0f, peptide.masses, peptide.intensities, emptyAAConstants);
			entries.add(entry);
		}
		Logger.logLine("Found "+entries.size()+" total peptide entries");

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
		
		LibraryFile library=new LibraryFile();
		library.openFile();
		Logger.logLine("Writing library file "+library.getName());
		library.dropIndices();
		library.addEntries(entries);
		library.addProteinsFromEntries(entries);
		library.createIndices();
		library.saveAsFile(libraryFile);
		return library;
	}
	
	public static IonType getIonType(String s) {
		if (s.length()==1) {
			switch (s.charAt(0)) {
			case 'A': return IonType.a;
			case 'B': return IonType.b;
			case 'C': return IonType.c;
			case 'X': return IonType.x;
			case 'Y': return IonType.y;
			case 'Z': return IonType.z;
			}
		} else if (s.charAt(1)=='2') {
			switch (s.charAt(0)) {
			case 'A': return IonType.ap2;
			case 'B': return IonType.bp2;
			case 'C': return IonType.cp2;
			case 'X': return IonType.xp2;
			case 'Y': return IonType.yp2;
			case 'Z': return IonType.zp2;
			}
		}
		
		Logger.errorLine("Found unexpected ion annotation ["+s+"]");
		throw new EncyclopediaException("Unexpected ion annotation ["+s+"]");
	}
	
	public static HashMap<String, SimplePeptidePrecursor> readInputFile(File peprecFile) {
		HashMap<String, SimplePeptidePrecursor> peptides=new HashMap<>();
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String id=row.get("spec_id");
				String ptmCodes=row.get("modifications");
				String peptideSequence=row.get("peptide");
				byte precursorCharge=Byte.parseByte(row.get("charge"));

				
				
				SimplePeptidePrecursor peptide=new SimplePeptidePrecursor(MS2PIPWriter.getPeptideModSeq(peptideSequence, ptmCodes), precursorCharge, emptyAAConstants);
				peptides.put(id, peptide);
			}

			@Override
			public void cleanup() {
			}
		};

		TableParser.parseSSV(peprecFile, muscle);
		return peptides;
	}
}
