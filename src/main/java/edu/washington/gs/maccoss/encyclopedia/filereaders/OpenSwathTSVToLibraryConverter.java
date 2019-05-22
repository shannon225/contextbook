package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideAccessionMatchingTrie;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class OpenSwathTSVToLibraryConverter {


	public static LibraryFile convertOpenSwathTSV(File tsvFile, File fastaFile, SearchParameters parameters) {
		String absolutePath=tsvFile.getAbsolutePath();
		File libraryFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		return convertFromOpenSwathTSV(tsvFile, fastaFile, libraryFile, parameters);
	}
	
	public static String getFromMap(Map<String, String> row, String... options) {
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
	
	public static void convertToOpenSwathTSV(SearchParameters params, final File elibFile, File tsvFile) throws IOException, SQLException, DataFormatException {
		LibraryFile library=new LibraryFile();
		library.openFile(elibFile);
		final AminoAcidConstants aaConstants=new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
		final ArrayList<LibraryEntry> allEntries=library.getAllEntries(false,  aaConstants);
		Logger.logLine("Found "+allEntries.size()+" entries from "+elibFile.getName()+". Writing to ["+tsvFile.getAbsolutePath()+"]...");

		try {
			PrintWriter writer=new PrintWriter(tsvFile, "UTF-8");
			HashSet<String> alreadyUsed=new HashSet<>();
			writer.println(General.toString(new String[] {
					"transition_group_id", "transition_name", "ProteinId", "PeptideSequence", "ModifiedPeptideSequence", "FullPeptideName", "RetentionTime", "PrecursorMz", "PrecursorCharge",
					"ProductMz", "ProductCharge", "LibraryIntensity", "FragmentIonType", "FragmentSeriesNumber", "IsDecoy", "quantifying_transition"}, "\t"));
			for (LibraryEntry e : allEntries) {
				AnnotatedLibraryEntry entry=AnnotatedLibraryEntry.getAnnotationsOnly(e, params);
				
				Object ModifiedSequence=PeptideUtils.formatForTPP(entry.getPeptideModSeq());
				Object transition_group_id=ModifiedSequence+"+"+entry.getPrecursorCharge();
				Object PeptideSequence=entry.getPeptideSeq();
				Object RetentionTime=new Float(entry.getRetentionTime());
				Object PrecursorMz=new Double(entry.getPrecursorMZ());
				Object PrecursorCharge=new Integer(entry.getPrecursorCharge());
				Object IsDecoy=entry.isDecoy()?"1":"0";
				Object ProteinId=PSMData.accessionsToString(entry.getAccessions());
				
				double[] masses=entry.getMassArray();
				float[] intensities=entry.getIntensityArray();
				FragmentIon[] ions=entry.getIonAnnotations();
				
				for (int i=0; i<ions.length; i++) {
					if (ions[i]!=null) {
						Object ProductCharge=new Byte(IonType.getCharge(ions[i].type));
						Object FragmentIonType=IonType.getType(ions[i].type);
						Object FragmentIonOrdinal=ions[i].index;
						Object transition_name=transition_group_id+"_"+FragmentIonType+FragmentIonOrdinal+"+"+PrecursorCharge;
						Object ProductMz=new Double(masses[i]);
						Object LibraryIntensity=new Float(intensities[i]);
						Object quantifying_transition="1";
						
						if (!alreadyUsed.contains(transition_name)) {
							alreadyUsed.add((String)transition_name);
							writer.println(General.toString(new Object[] {
									transition_group_id,transition_name,ProteinId,PeptideSequence,ModifiedSequence,ModifiedSequence,RetentionTime,PrecursorMz,PrecursorCharge,
									ProductMz,ProductCharge,LibraryIntensity,FragmentIonType,FragmentIonOrdinal,IsDecoy, quantifying_transition}, "\t"));
						}
					}
				}
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			Logger.errorLine("Error writing library to OpenSWATH file.");
			Logger.errorException(e);
		}
		
		library.close();
		Logger.logLine("Finished reading "+tsvFile.getName());
	}

	public static LibraryInterface convertFromOpenSwathTSV(File tsvFile, File fastaFile, SearchParameters parameters) {
		String absolutePath=tsvFile.getAbsolutePath();
		File libraryFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		return convertFromOpenSwathTSV(tsvFile, fastaFile, libraryFile, parameters);
	}

	public static LibraryFile convertFromOpenSwathTSV(File tsvFile, File fastaFile, File libraryFile, SearchParameters parameters) {
		AminoAcidConstants aaConstants=parameters.getAAConstants();
		try {
			final ArrayList<ImmutablePeptideEntry> peptides=new ArrayList<ImmutablePeptideEntry>();
			TableParserMuscle muscle=new TableParserMuscle() {
				private PeptideEntry lastPeptide=null;
				private String lastGroup=null;
				@Override
				public void processRow(Map<String, String> row) {
					int decoy=Integer.parseInt(row.get("decoy"));
					if (decoy!=0) return;
					
					String group=row.get("transition_group_id");
					String peptideModSeq=parseMods(getFromMap(row, "ModifiedPeptideSequence", "FullUniModPeptideName", "FullPeptideName", "ModifiedSequence", "PeptideSequence", "Sequence", "StrippedSequence"));
					byte charge=Byte.parseByte(row.get("PrecursorCharge"));
					double productMz=Double.parseDouble(getFromMap(row, "ProductMz", "FragmentMz"));
					float libraryIntensity=Float.parseFloat(getFromMap(row, "LibraryIntensity", "RelativeFragmentIntensity"));
					float iRT=Float.parseFloat(getFromMap(row, "NormalizedRetentionTime", "RetentionTime", "Tr_recalibrated", "iRT", "RetentionTimeCalculatorScore"));
					
					
					if (!group.equals(lastGroup)) {
						if (lastPeptide!=null) peptides.add(new ImmutablePeptideEntry(lastPeptide));
						lastGroup=group;
						lastPeptide=new PeptideEntry(peptideModSeq, charge, iRT);
					}
					lastPeptide.addPeak(new Peak(productMz, libraryIntensity));
				}
				
				@Override
				public void cleanup() {
					if (lastPeptide!=null) peptides.add(new ImmutablePeptideEntry(lastPeptide));
				}
			};
			
			TableParser.parseTSV(tsvFile, muscle);

			return processPeptideEntries(tsvFile.getName(), fastaFile, libraryFile, parameters, aaConstants, peptides);

		} catch (Exception e) {
			Logger.errorLine("Error parsing OpenSwath TSV:");
			Logger.errorException(e);
			throw new EncyclopediaException(e);
		}
	}

	public static LibraryFile processPeptideEntries(String sourceFile, File fastaFile, File libraryFile, SearchParameters parameters, AminoAcidConstants aaConstants,
			final ArrayList<ImmutablePeptideEntry> peptides) throws IOException, SQLException {
		ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
		for (ImmutablePeptideEntry peptide : peptides) {
			double precursorMZ=aaConstants.getChargedMass(peptide.peptideModSeq, peptide.charge);
			HashSet<String> accessions=new HashSet<>();
			
			if (fastaFile==null) {
				accessions.add(PeptideUtils.getPeptideSeq(peptide.peptideModSeq));
			}
			
			LibraryEntry entry=new LibraryEntry(sourceFile, accessions, precursorMZ, peptide.charge, peptide.peptideModSeq, 1, peptide.rt, 0.0f, peptide.masses, peptide.intensities, aaConstants);
			entries.add(entry);
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

}
