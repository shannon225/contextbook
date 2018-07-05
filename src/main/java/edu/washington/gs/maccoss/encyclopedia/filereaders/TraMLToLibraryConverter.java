package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.hupo.psi.ms.traml.CvParamType;
import org.hupo.psi.ms.traml.ModificationType;
import org.hupo.psi.ms.traml.PeptideType;
import org.hupo.psi.ms.traml.RetentionTimeType;
import org.hupo.psi.ms.traml.TraMLType;
import org.hupo.psi.ms.traml.TransitionType;
import org.systemsbiology.apps.tramlparser.TraMLParser;

import com.compomics.jtraml.enumeration.FrequentOBoEnum;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideAccessionMatchingTrie;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;

public class TraMLToLibraryConverter {
	private static org.apache.log4j.Logger logger=org.apache.log4j.Logger.getLogger(TraMLToLibraryConverter.class);

	public static void main(String[] args) {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		File tramlFile=new File("/Users/searleb/Documents/phospho_localization/rosenberger/dda_psgs_consensus_spectral_library/psgs_standard_consensus_filtered.TraML");
		//File fastaFile=new File("/Users/searleb/Documents/phospho_localization/rosenberger/PSGS_reference_peptides_decoys.fasta");
		Long time=System.currentTimeMillis();
		convertTraML(tramlFile, null, params);
		System.out.println(System.currentTimeMillis()-time);
	}

	public static LibraryFile convertTraML(File tramlFile, File fastaFile, SearchParameters parameters) {
		String absolutePath=tramlFile.getAbsolutePath();
		File libraryFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		return convertTraML(tramlFile, fastaFile, libraryFile, parameters);
	}

	public static LibraryFile convertTraML(File tramlFile, File fastaFile, File libraryFile, SearchParameters parameters) {
		AminoAcidConstants aaConstants=parameters.getAAConstants();
		TraMLParser traMLParser=new TraMLParser();

		HashMap<String, PeptideEntry> peptideByID=new HashMap<>();

		try {
			traMLParser.parse_file(tramlFile.getCanonicalPath(), logger);
			TraMLType traML=traMLParser.getTraML();
			List<TransitionType> transitionTypeList=traML.getTransitionList().getTransition();

			for (TransitionType transitionType : transitionTypeList) {
				PeptideType peptideType=(PeptideType)transitionType.getPeptideRef();

				PeptideEntry peptide=peptideByID.get(peptideType.getId());
				if (peptide==null) {
					String peptideModSeq=getPeptideModSeq(peptideType);
					byte charge=0;
					for (CvParamType param : peptideType.getCvParam()) {
						if ("charge state".equals(param.getName())) {
							charge=Byte.parseByte(param.getValue());
							break;
						}
					}
					float rt=getRT(peptideType);
					peptide=new PeptideEntry(peptideModSeq, charge, rt);
					peptideByID.put(peptideType.getId(), peptide);
				}

				float intensity=0.0f;
				double mass=0.0;
				for (CvParamType param : transitionType.getCvParam()) {
					if ("product ion intensity".equals(param.getName())) {
						intensity=Float.parseFloat(param.getValue());
					}
				}
				for (CvParamType param : transitionType.getProduct().getCvParam()) {
					if ("isolation window target m/z".equals(param.getName())) {
						mass=Double.parseDouble(param.getValue());
					}
				}
				peptide.addPeak(new Peak(mass, intensity));

			}

			ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
			for (PeptideEntry peptide : peptideByID.values()) {
				Collections.sort(peptide.peaks);
				Pair<double[], float[]> peakArrays=Peak.toArrays(peptide.peaks);
				double precursorMZ=aaConstants.getChargedMass(peptide.peptideModSeq, peptide.charge);
				HashSet<String> accessions=new HashSet<>();
				
				if (fastaFile==null) {
					accessions.add(PeptideUtils.getPeptideSeq(peptide.peptideModSeq));
				}
				
				LibraryEntry entry=new LibraryEntry(tramlFile.getName(), accessions, precursorMZ, peptide.charge, peptide.peptideModSeq, 1, peptide.rt, 0.0f, peakArrays.x, peakArrays.y, aaConstants);
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

		} catch (Exception e) {
			Logger.errorLine("Error parsing TraML:");
			Logger.errorException(e);
			throw new EncyclopediaException(e);
		}
	}

	private static float getRT(PeptideType peptideType) {
		for (RetentionTimeType rtType : peptideType.getRetentionTimeList().getRetentionTime()) {
			List<CvParamType> params=rtType.getCvParam();
			for (CvParamType param : params) {
				if (FrequentOBoEnum.RETENTION_TIME.getName().equals(param.getName())) {
					return  Float.parseFloat(param.getValue())*60f;
				} else if (FrequentOBoEnum.RETENTION_TIME_NORMALIZED.getName().equals(param.getName())) {
					return Float.parseFloat(param.getValue())*60f;
				}
			}
		}
		return 0.0f;
	}

	private static String getPeptideModSeq(PeptideType peptideType) {
		try {
			String sequence=peptideType.getSequence();
			String[] aas=new String[sequence.length()];
			for (int i=0; i<aas.length; i++) {
				aas[i]=Character.toString(sequence.charAt(i));
			}
			List<ModificationType> modTypes=peptideType.getModification();

			for (ModificationType modType : modTypes) {
				Double modMass=modType.getMonoisotopicMassDelta();
				int index=modType.getLocation()-1;
				if (index==aas.length) index=aas.length-1;
				aas[index]=aas[index]+"["+(modMass>=0?"+":"")+modMass+"]";
			}

			StringBuilder sb=new StringBuilder();
			for (int i=0; i<aas.length; i++) {
				sb.append(aas[i]);
			}
			String peptideModSeq=sb.toString();
			return peptideModSeq;
		} catch (Exception e) {
			Logger.errorLine("TraML parsing error! Trying to parse "+peptideType.getSequence()+", ["+peptideType.getId()+"]. "+peptideType.getModification().size()+" mods:");
			for (ModificationType modType : peptideType.getModification()) {
				Logger.errorLine("\t"+modType.getAverageMassDelta()+" (index: "+(modType.getLocation()-1)+")");
			}
			throw new EncyclopediaException(e);
		}
	}

}
