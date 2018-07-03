package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Map;
import java.util.StringTokenizer;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideAccessionMatchingTrie;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;

public class OpenSwathTSVToLibraryConverter {


	public static LibraryFile convertOpenSwathTSV(File tsvFile, File fastaFile, AminoAcidConstants aaConstants) {
		String absolutePath=tsvFile.getAbsolutePath();
		File libraryFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		return convertOpenSwathTSV(tsvFile, fastaFile, libraryFile, aaConstants);
	}
	
	private static String getFromMap(Map<String, String> row, String... options) {
		for (String option : options) {
			String value=row.get(option);
			if (value!=null) return value;
		}
		return null;
	}
	
	private static String parseMods(String structuredSequence) {
		
		if (structuredSequence.indexOf('(')>=0) {
			 // Unimod: .(UniMod:1)PEPC(UniMod:4)PEPM(UniMod:35)PEPR.(UniMod:2)
			// FIXME Parsing Unimod inside of OpenSwath TSVs isn't supported yet!
			throw new EncyclopediaException("Parsing Unimod inside of OpenSwath TSVs isn't supported yet!");
		} else {
			if (structuredSequence.indexOf('[')>=0) {
				 // TPP:    n[43]PEPC[160]PEPM[147]PEPRc[16]
				StringBuilder sb=new StringBuilder(structuredSequence);
				final int nIndex=structuredSequence.indexOf('n');
				if (nIndex>=0) sb=sb.deleteCharAt(nIndex);
				final int cIndex=structuredSequence.indexOf('c');
				if (cIndex>=0) sb=sb.deleteCharAt(cIndex);
				return sb.toString();
				
			} else if (structuredSequence.indexOf('.')>=0) {
				 // Unimod: .(UniMod:1)PEPC(UniMod:4)PEPM(UniMod:35)PEPR.(UniMod:2) (but no mods)
				StringTokenizer st=new StringTokenizer(structuredSequence, ".");
				st.nextToken();
				return st.nextToken();
				
			} else if (structuredSequence.indexOf('n')>=0) {
				 // TPP:    n[43]PEPC[160]PEPM[147]PEPRc[16] (but no mods)
				return structuredSequence.replace('n', ' ').replace('c', ' ').trim();
				
			} else {
				return structuredSequence;
			}
		}
	}

	public static LibraryFile convertOpenSwathTSV(File tsvFile, File fastaFile, File libraryFile, AminoAcidConstants aaConstants) {
		try {
			final ArrayList<PeptideEntry> peptides=new ArrayList<PeptideEntry>();
			TableParserMuscle muscle=new TableParserMuscle() {
				private PeptideEntry lastPeptide=null;
				private String lastGroup=null;
				@Override
				public void processRow(Map<String, String> row) {
					int decoy=Integer.parseInt("decoy");
					if (decoy!=0) return;
					
					String group=row.get("transition_group_id");
					String peptideModSeq=parseMods(getFromMap(row, "ModifiedPeptideSequence", "FullUniModPeptideName", "FullPeptideName", "ModifiedSequence", "PeptideSequence", "Sequence", "StrippedSequence"));
					byte charge=Byte.parseByte("PrecursorCharge");
					double productMz=Double.parseDouble(getFromMap(row, "ProductMz", "FragmentMz"));
					float libraryIntensity=Float.parseFloat(getFromMap(row, "LibraryIntensity", "RelativeFragmentIntensity"));
					float iRT=Float.parseFloat(getFromMap(row, "NormalizedRetentionTime", "RetentionTime", "Tr_recalibrated", "iRT", "RetentionTimeCalculatorScore"));
					
					
					if (lastGroup!=group) {
						if (lastPeptide!=null) peptides.add(lastPeptide);
						
						lastPeptide=new PeptideEntry(peptideModSeq, charge, iRT);
					}
					lastPeptide.addPeak(new Peak(productMz, libraryIntensity));
				}
				
				@Override
				public void cleanup() {
					if (lastPeptide!=null) peptides.add(lastPeptide);
				}
			};
			
			TableParser.parseTSV(tsvFile, muscle);

			ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
			for (PeptideEntry peptide : peptides) {
				Collections.sort(peptide.peaks);
				Pair<double[], float[]> peakArrays=Peak.toArrays(peptide.peaks);
				double precursorMZ=aaConstants.getChargedMass(peptide.peptideModSeq, peptide.charge);
				HashSet<String> accessions=new HashSet<>();
				
				if (fastaFile==null) {
					accessions.add(PeptideUtils.getPeptideSeq(peptide.peptideModSeq));
				}
				
				LibraryEntry entry=new LibraryEntry(tsvFile.getName(), accessions, precursorMZ, peptide.charge, peptide.peptideModSeq, 1, peptide.rt, 0.0f, peakArrays.x, peakArrays.y, aaConstants);
				entries.add(entry);
			}

			if (fastaFile!=null) {
				Logger.logLine("Reading Fasta file "+fastaFile.getName());
				ArrayList<FastaEntryInterface> proteins=FastaReader.readFasta(fastaFile);
			
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

		} catch (Exception e) {
			Logger.errorLine("Error parsing OpenSwath TSV:");
			Logger.errorException(e);
			throw new EncyclopediaException(e);
		}
	}

}
