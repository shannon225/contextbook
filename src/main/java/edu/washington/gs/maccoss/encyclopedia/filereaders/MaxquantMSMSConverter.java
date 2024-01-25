package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Optional;
import java.util.StringTokenizer;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideAccessionMatchingTrie;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.procedure.TCharDoubleProcedure;
import gnu.trove.procedure.TObjectDoubleProcedure;

public class MaxquantMSMSConverter {
	private static final String VALUE_LIST_DELIMITER = ";";
	
	public static void main(String[] args) {
		File tsvFile=new File("/Users/searleb/Documents/swaney/CID_vs_HCD_enzymes/CombinedLibrarySearch_MQ/CID/txt/msms.txt");
		File fastaFile=new File("/Users/searleb/Documents/swaney/CID_vs_HCD_enzymes/Uniprot_human_9606_canonical_020320.fasta");
		File libraryFile=new File("/Users/searleb/Documents/swaney/CID_vs_HCD_enzymes/CombinedLibrarySearch_MQ/CID/txt/msms.dlib");
		convertFromMSMSTSV(tsvFile, fastaFile, libraryFile, SearchParameterParser.getDefaultParametersObject());
	}
	
	/**
	 * FIXME FRAGILE! FIGURE OUT HOW THEY REALLY ANNOTATE MODS BEYOND THESE
	 */
	public static final TObjectDoubleHashMap<String> modMap=new TObjectDoubleHashMap<>();
	public static void initPTMs() {
		modMap.clear();
		modMap.put("ac", 42.010565);
		modMap.put("cam", 57.0214635);
		modMap.put("p", 79.966331);
		modMap.put("ph", 79.966331);
		modMap.put("me", 14.015650);
		modMap.put("ub", 114.042927);
		modMap.put("ox", 15.994915);
		modMap.put("Phospho (STY)", 79.966331);
		modMap.put("Oxidation (M)", 15.994915);
		modMap.put("Acetyl (Protein N-term)", 42.010565);
		modMap.put("Acetyl (N-term)", 42.010565);
		modMap.put("Acetyl (K)", 42.010565);
		modMap.put("Carbamidomethyl (C)", 57.0214635);
	}
	
	private static TCharDoubleHashMap checkForFixedModsFromMSMS(File tsvFile, SearchParameters parameters) {
		File summaryFile=new File(tsvFile.getParentFile(), "summary.txt");
		if (summaryFile.exists()&&summaryFile.isFile()&&summaryFile.canRead()) {
			TCharDoubleHashMap fixedMods=new TCharDoubleHashMap();
			TableParserMuscle muscle=new TableParserMuscle() {
				@Override
				public void processRow(Map<String, String> row) {
					String fixedModText=row.get("Fixed modifications");
					if (fixedModText!=null) {
						StringTokenizer st=new StringTokenizer(fixedModText, ";");
						while (st.hasMoreTokens()) {
							String name=st.nextToken();
							if (modMap.contains(name)) {
								double mass=modMap.get(name);
								
								int startIndex=name.lastIndexOf('(')+1;
								String aas=name.substring(startIndex, name.lastIndexOf(')'));
								
								for (char aa : aas.toCharArray()) {
									fixedMods.put(aa, mass);
								}
							}
						}
					}
					
					// we don't do anything with this section of the file -- it's just to add log messages for expected errors
					String variableModText=row.get("Variable modifications");
					if (variableModText!=null) {
						StringTokenizer st=new StringTokenizer(variableModText, ";");
						while (st.hasMoreTokens()) {
							String name=st.nextToken();
							if (!modMap.contains(name)) {
								Logger.errorLine("Unexpected PTM from MaxQuant: ["+name+"], sorry peptides modified in this way will not be parsed correctly!");
							}
						}
					}
				}
				
				@Override
				public void cleanup() {
				}
			};
			TableParser.parseTSV(summaryFile, muscle);
			
			return fixedMods;
		} else {
			Logger.errorLine("Could not find MaxQuant summary file in "+tsvFile.getParentFile().getAbsolutePath()+", proceeding without fixed PTMs.");
			return parameters.getAAConstants().getFixedMods();
		}
	}

	public static LibraryFile convertFromMSMSTSV(File tsvFile, File fastaFile, File libraryFile, SearchParameters parameters) {
		initPTMs();
		TCharDoubleHashMap fixedMods=checkForFixedModsFromMSMS(tsvFile, parameters);
		
		AminoAcidConstants aaConstants=parameters.getAAConstants();
		try {
			final HashMap<String, ScoredObject<ImmutablePeptideEntry>> peptides=new HashMap<>();
			TableParserMuscle muscle=new TableParserMuscle() {
				@Override
				public void processRow(Map<String, String> row) {
					String rawFile=row.get("Raw file");
					if (rawFile==null) rawFile=tsvFile.getName();
					String peptideModSeq=parseMods(row.get("Modified sequence"), fixedMods);
					byte charge=Byte.parseByte(row.get("Charge"));
					float[] intensities=parseFloats(row.get("Intensities"));
					double[] masses=parseDoubles(row.get("Masses"));

					// Compatible with new versions: MaxQuant 1.6.10+
					// Note that the name changed between 1.6.3 and 1.6.6
					double[] massDeviations=parseDoubles(row.containsKey("Mass Deviations [Da]") ? row.get("Mass Deviations [Da]") : row.get("Mass deviations [Da]"));

					float rt=Float.parseFloat(row.get("Retention time"));

					// Compatible with new versions: MaxQuant 1.6.10+
					// Note that the name changed between 1.6.3 and 1.6.6
					float rtOffset=Float.parseFloat(row.containsKey("Precursor Apex Offset Time") ? row.get("Precursor Apex Offset Time") : row.get("Precursor apex offset time"));

					float score=(float)-Log.protectedLog10(Double.parseDouble(row.get("PEP")));

					if (Float.isNaN(rtOffset)) {
						rtOffset=0.0f;
					}
					if ((!Float.isFinite(rt))||(!Float.isFinite(rtOffset))) {
						System.err.println("RT ERROR: rt:"+rt+", rtOffset:"+rtOffset);
					}

					// RTs are in mins
					PeptideEntry pep=new PeptideEntry(peptideModSeq, charge, (rt+rtOffset)*60f, rawFile);
					for (int i = 0; i < masses.length; i++) {
						pep.addPeak(new Peak(masses[i]+massDeviations[i], intensities[i]));
					}
					String key = peptideModSeq+"+"+charge;
					ScoredObject<ImmutablePeptideEntry> value = peptides.get(key);
					if (value==null||value.getScore()<score) {
						value = new ScoredObject<ImmutablePeptideEntry>(score, new ImmutablePeptideEntry(pep));
						peptides.put(key, value);
					}
				}
				
				@Override
				public void cleanup() {
				}
			};
			
			TableParser.parseTSV(tsvFile, muscle);
			
			return processPeptideEntries(tsvFile.getName(), fastaFile, libraryFile, parameters, aaConstants, peptides.values());

		} catch (Exception e) {
			Logger.errorLine("Error parsing Maxquant msms.txt:");
			Logger.errorException(e);
			throw e instanceof EncyclopediaException
					? (EncyclopediaException) e
					: new EncyclopediaException(e);
		}
	}
	
	private static class StringWrapper {
		String s;
		public StringWrapper(String s) {
			this.s = s;
		}
	}
	
	/**
	 * @param seq
	 * @return
	 */
	public static String parseMods(String seq, TCharDoubleHashMap fixedMods) {
		seq=seq.replace("_", "");

		// forEachEntry is threadsafe, so this is dumb, but necessary to modify the "seq" object
		final StringWrapper seqWrapper=new StringWrapper(seq);
		
		modMap.forEachEntry(new TObjectDoubleProcedure<String>() {
			@Override
			public boolean execute(String a, double b) {
				seqWrapper.s=seqWrapper.s.replace("("+a+")", "["+b+"]");
				return true;
			}
		});
		
		// fixed PTMs get added second, after the longer strings get adjusted
		fixedMods.forEachEntry(new TCharDoubleProcedure() {
			@Override
			public boolean execute(char a, double b) {
				if (a==AminoAcidConstants.N_TERM) {
					seqWrapper.s="["+b+"]"+seqWrapper.s;
				} else if (a==AminoAcidConstants.C_TERM) {
					seqWrapper.s=seqWrapper.s+"["+b+"]";
				} else {
					seqWrapper.s=seqWrapper.s.replace(a+"", a+"["+b+"]");
				}
				return true;
			}
		});
		return seqWrapper.s;
	}

	public static LibraryFile processPeptideEntries(String sourceFile, File fastaFile, File libraryFile, SearchParameters parameters, AminoAcidConstants aaConstants,
			final Collection<ScoredObject<ImmutablePeptideEntry>> peptides) throws IOException, SQLException {
		ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
		for (ScoredObject<ImmutablePeptideEntry> scoredObj : peptides) {
			ImmutablePeptideEntry peptide=scoredObj.y;
			double precursorMZ=aaConstants.getChargedMass(peptide.peptideModSeq, peptide.charge);
			HashSet<String> accessions=new HashSet<>();
			
			if (fastaFile==null) {
				accessions.add(PeptideUtils.getPeptideSeq(peptide.peptideModSeq));
			}
			
			LibraryEntry entry=new LibraryEntry(peptide.sourceFile, accessions, precursorMZ, peptide.charge, peptide.peptideModSeq, 1, peptide.rt, (float)Math.pow(10, -scoredObj.x), peptide.masses, peptide.intensities, Optional.empty(), // FIXME add ion mobility to parser
					aaConstants);
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

	private static float[] parseFloats(String values) {
		if (null == values || values.isEmpty()) {
			return new float[0];
		} else if (!values.contains(VALUE_LIST_DELIMITER)) {
			return parseFloats(new String[]{values});
		} else {
			return parseFloats(values.split(VALUE_LIST_DELIMITER));
		}
	}

	private static float[] parseFloats(String[] values) {
		float[] f=new float[values.length];
		for (int i = 0; i < f.length; i++) {
			f[i]=Float.parseFloat(values[i]);
		}
		return f;
	}

	private static double[] parseDoubles(String values) {
		if (null == values || values.isEmpty()) {
			return new double[0];
		} else if (!values.contains(VALUE_LIST_DELIMITER)) {
			return parseDoubles(new String[]{values});
		} else {
			return parseDoubles(values.split(VALUE_LIST_DELIMITER));
		}
	}

	private static double[] parseDoubles(String[] values) {
		double[] f=new double[values.length];
		for (int i = 0; i < f.length; i++) {
			f[i]=Double.parseDouble(values[i]);
		}
		return f;
	}

}
