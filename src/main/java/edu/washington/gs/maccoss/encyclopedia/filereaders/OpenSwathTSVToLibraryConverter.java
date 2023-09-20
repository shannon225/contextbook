package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.zip.DataFormatException;

import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideAccessionMatchingTrie;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PTMMap.PostTranslationalModification;
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
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.procedure.TObjectProcedure;

public class OpenSwathTSVToLibraryConverter {
	public static final List<String> PEPTIDE_HEADERS = ImmutableList.copyOf(new String[] {
			"ModifiedSequence",
			"FullPeptideName",
			"FullUniModPeptideName",
			"LabeledSequence",
			"ModifiedPeptide",
			"ModifiedPeptideSequence",
			"PeptideSequence",
			"Sequence",
			"StrippedSequence",
			"modification_sequence",
	});
	public static final List<String> CHARGE_HEADERS = ImmutableList.copyOf(new String[] {
			"PrecursorCharge",
			"Charge",
			"prec_z",
	});
	public static final List<String> FRAG_MZ_HEADERS = ImmutableList.copyOf(new String[] {
			"ProductMz",
			"FragmentMz",
	});
	public static final List<String> FRAG_INTEN_HEADERS = ImmutableList.copyOf(new String[] {
			"LibraryIntensity",
			"RelativeIntensity",
			"RelativeFragmentIntensity",
			"RelativeFragmentIonIntensity",
			"relative_intensity",
	});
	public static final List<String> RT_HEADERS = ImmutableList.copyOf(new String[] {
			"NormalizedRetentionTime",
			"RetentionTime",
			"RetentionTimeCalculatorScore",
			"RT_detected",
			"Tr_recalibrated",
			"iRT",
	});
	public static final List<String> TRANSITION_GROUP_HEADERS = ImmutableList.copyOf(new String[] {
			"TransitionGroupId",
			"ElutionGroup",
			"transition_group_id",
	});

	public static LibraryFile convertOpenSwathTSV(File tsvFile, File fastaFile, SearchParameters parameters) {
		String absolutePath=tsvFile.getAbsolutePath();
		File libraryFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		return convertFromOpenSwathTSV(tsvFile, fastaFile, libraryFile, parameters);
	}

	public static String getFromMap(Map<String, String> row, Collection<? extends String> options) {
		return getFromMap(row, options.toArray(new String[0]));
	}
	
	public static String getFromMap(Map<String, String> row, String... options) {
		for (String option : options) {
			String value=row.get(option);
			if (value!=null) return value;
		}
		return null;
	}
	
	static String parseMods(String structuredSequence) {
		if (structuredSequence.indexOf('(')>=0) {
			structuredSequence=structuredSequence.replace(".", "");
			// Unimod: .(UniMod:1)PEPC(UniMod:4)PEPM(UniMod:35)PEPR.(UniMod:2)

			char[] ca=structuredSequence.toCharArray();
			
			ArrayList<String> aas=new ArrayList<String>();
			for (int i = 0; i < ca.length; i++) {
				if (ca[i]=='(') {
					StringBuilder sb=new StringBuilder();
					i++;
					while (ca[i]!=')') {
						sb.append(ca[i]);
						i++;
					}
					if (aas.size()==0) {
						// handling of n-termini mods assumes you can't have multiple []s in a row
						i++;
						aas.add(Character.toString(ca[i]));
					}
					
					PostTranslationalModification ptm=PTMMap.getPTM(sb.toString().toUpperCase());
					double modificationMass = ptm.getDeltaMass();
					
					String aaString=aas.get(aas.size()-1);
					aas.set(aas.size()-1, aaString+(modificationMass>=0?"[+":"[")+modificationMass+"]");
				} else {
					aas.add(Character.toString(ca[i]));
				}
			}
			
			StringBuilder sb=new StringBuilder();
			for (String aa : aas) {
				sb.append(aa);
			}
			return sb.toString();
			
		} else {
			if (structuredSequence.indexOf('[')>=0) {
				 // TPP:    n[43]PEPC[160]PEPM[147]PEPRc[16]
				StringBuilder sb=new StringBuilder(structuredSequence);
				final int nIndex=structuredSequence.indexOf(AminoAcidConstants.N_TERM);
				if (nIndex>=0) sb=sb.deleteCharAt(nIndex);
				final int cIndex=sb.toString().indexOf(AminoAcidConstants.C_TERM); // need to get index in possibly-modified string!
				if (cIndex>=0) sb=sb.deleteCharAt(cIndex);
				return sb.toString();
				
			} else if (structuredSequence.indexOf('.')>=0) {
				 // Unimod: .(UniMod:1)PEPC(UniMod:4)PEPM(UniMod:35)PEPR.(UniMod:2) (but no mods)
				StringTokenizer st=new StringTokenizer(structuredSequence, ".");
				return st.nextToken();
				
			} else if (structuredSequence.indexOf(AminoAcidConstants.N_TERM)>=0) {
				 // TPP:    n[43]PEPC[160]PEPM[147]PEPRc[16] (but no mods)
				return structuredSequence.replace(AminoAcidConstants.N_TERM, ' ').replace(AminoAcidConstants.C_TERM, ' ').trim();
				
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
						Object ProductCharge=new Byte(IonType.getCharge(ions[i].getType()));
						Object FragmentIonType=IonType.getType(ions[i].getType());
						Object FragmentIonOrdinal=ions[i].getIndex();
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

	/**
	 * @return A <emph>closed</emph> {@code LibraryFile} instance pointing to {@code elibFile}
	 *         containing the results of conversion.
	 */
	public static LibraryFile convertFromOpenSwathTSV(File tsvFile, File fastaFile, SearchParameters parameters) {
		String absolutePath=tsvFile.getAbsolutePath();
		File libraryFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		return convertFromOpenSwathTSV(tsvFile, fastaFile, libraryFile, parameters);
	}

	/**
	 * @return A <emph>closed</emph> {@code LibraryFile} instance pointing to {@code elibFile}
	 *         containing the results of conversion.
	 */
	public static LibraryFile convertFromOpenSwathTSV(File tsvFile, File fastaFile, File libraryFile, SearchParameters parameters) {
		String sourceFile = tsvFile.getName();
		AminoAcidConstants aaConstants=parameters.getAAConstants();
		try {
			final ArrayList<ImmutablePeptideEntry> peptides=new ArrayList<>();
			final TIntObjectHashMap<PeptideEntry> peptideMap=new TIntObjectHashMap<>();
			TableParserMuscle muscle=new TableParserMuscle() {
				/**
				 * Used only if the {@code transition_group_id} column is missing.
				 */
				int peptideCount = 0;

				/**
				 * Used only if the {@code transition_group_id} column is missing.
				 */
				String lastPeptideModSeq = null;

				/**
				 * Used only if the {@code transition_group_id} column is missing.
				 */
				byte lastPeptideCharge = -1;

				@Override
				public void processRow(Map<String, String> row) {
					String decoy=getFromMap(row, "decoy", "Decoy");
					if (decoy!=null&&Integer.parseInt(decoy)!=0) return;

					String peptideModSeq=parseMods(getFromMap(row, PEPTIDE_HEADERS));
					byte charge=Byte.parseByte(getFromMap(row, CHARGE_HEADERS));
					double productMz=Double.parseDouble(getFromMap(row, FRAG_MZ_HEADERS));
					float libraryIntensity=Float.parseFloat(getFromMap(row, FRAG_INTEN_HEADERS));
					float iRT=Float.parseFloat(getFromMap(row, RT_HEADERS));

					final String groupIdString = getFromMap(row, TRANSITION_GROUP_HEADERS);

					final int group;
					if (null == groupIdString) {
						// Group IDs not reported; assume we can use peptideModSeq/charge. Requires that
						// each (peptideModSeq, charge) pair has only a single entry in the file and that
						// all the fragments for a given pair are reported on consecutive lines.
						if (null == lastPeptideModSeq || !lastPeptideModSeq.equals(peptideModSeq) || 0 > lastPeptideCharge || lastPeptideCharge != charge) {
							lastPeptideModSeq = peptideModSeq;
							lastPeptideCharge = charge;
							group = ++peptideCount;
						} else {
							group = peptideCount;
						}
					} else {
						group = Integer.parseInt(groupIdString);
					}

					PeptideEntry thisPeptide = peptideMap.get(group);
					if (thisPeptide==null) {
						thisPeptide=new PeptideEntry(peptideModSeq, charge, iRT*60f, sourceFile);
						peptideMap.put(group, thisPeptide);
					}

					if (libraryIntensity>0) {
						thisPeptide.addPeak(new Peak(productMz, libraryIntensity));
					}
				}
				
				@Override
				public void cleanup() {
					peptideMap.forEachValue(new TObjectProcedure<PeptideEntry>() {
						@Override
						public boolean execute(PeptideEntry pep) {
							peptides.add(new ImmutablePeptideEntry(pep));
							return true;
						}
					});
				}
			};
			
			TableParser.parseTSV(tsvFile, muscle);

			return processPeptideEntries(sourceFile, fastaFile, libraryFile, parameters, aaConstants, peptides);

		} catch (Exception e) {
			Logger.errorLine("Error parsing OpenSwath TSV:");
			Logger.errorException(e);
			throw new EncyclopediaException(e);
		}
	}

	/**
	 * @return A <emph>closed</emph> {@code LibraryFile} instance pointing to {@code elibFile}
	 *         containing the results of conversion.
	 */
	public static LibraryFile processPeptideEntries(String sourceFile, File fastaFile, File libraryFile, SearchParameters parameters, AminoAcidConstants aaConstants,
			final ArrayList<ImmutablePeptideEntry> peptides) throws IOException, SQLException {
		ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
		for (ImmutablePeptideEntry peptide : peptides) {
			double precursorMZ=aaConstants.getChargedMass(peptide.peptideModSeq, peptide.charge);
			HashSet<String> accessions=new HashSet<>();
			
			if (fastaFile==null) {
				accessions.add(PeptideUtils.getPeptideSeq(peptide.peptideModSeq));
			}
			
			LibraryEntry entry=new LibraryEntry(peptide.sourceFile, accessions, precursorMZ, peptide.charge, peptide.peptideModSeq, 1, peptide.rt, 0.0f, peptide.masses, peptide.intensities, aaConstants);
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
		library.close();
		return library;
	}

}
