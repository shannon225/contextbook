package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

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

public class MaxquantMSMSConverter {
	private static final String VALUE_LIST_DELIMITER = ";";

	public static LibraryFile convertFromMSMSTSV(File tsvFile, File fastaFile, File libraryFile, SearchParameters parameters) {
		AminoAcidConstants aaConstants=parameters.getAAConstants();
		try {
			final HashMap<String, ScoredObject<ImmutablePeptideEntry>> peptides=new HashMap<>();
			TableParserMuscle muscle=new TableParserMuscle() {
				@Override
				public void processRow(Map<String, String> row) {
					String peptideModSeq=parseMods(row.get("Modified sequence"));
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
					PeptideEntry pep=new PeptideEntry(peptideModSeq, charge, (rt+rtOffset)*60f);
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
			throw new EncyclopediaException(e);
		}
	}
	
	/**
	 * FIXME FRAGILE! FIGURE OUT HOW THEY REALLY ANNOTATE MODS BEYOND THESE
	 * @param seq
	 * @return
	 */
	public static String parseMods(String seq) {
		seq=seq.replace("_", "");
		seq=seq.replace("(ac)", "[42.010565]");
		seq=seq.replace("(cam)", "[57.0214635]");
		seq=seq.replace("(p)", "[79.966331]");
		seq=seq.replace("(ph)", "[79.966331]");
		seq=seq.replace("(me)", "[14.015650]");
		seq=seq.replace("(ub)", "[114.042927]");
		seq=seq.replace("(ox)", "[15.994915]");
		return seq;
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
			
			LibraryEntry entry=new LibraryEntry(sourceFile, accessions, precursorMZ, peptide.charge, peptide.peptideModSeq, 1, peptide.rt, (float)Math.pow(10, -scoredObj.x), peptide.masses, peptide.intensities, aaConstants);
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
