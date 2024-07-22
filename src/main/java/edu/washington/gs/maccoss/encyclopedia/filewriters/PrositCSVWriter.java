package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.zip.DataFormatException;

import org.apache.commons.lang3.StringUtils;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SimplePeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class PrositCSVWriter {
	public static void main4(String[] args) throws Exception {
		File fileName=new File("/Users/searleb/Documents/grants/phosphotau mutagenesis/barcodes09_prosit.csv");
		PrintWriter writer=new PrintWriter(fileName);
		// removing R, K, M, C, G, F, Y, W, and I
		char[] aas="ADEHLNQTV".toCharArray();
		//char[] aas="ADEHLNPQSTVRKMCGFYW".toCharArray();
		byte nce=33;
		byte charge=2;
		writer.println("modified_sequence,collision_energy,precursor_charge");
		
		for (int i = 0; i < aas.length; i++) {
			for (int j = 0; j < aas.length; j++) {
				for (int k = 0; k < aas.length; k++) {
					for (int l = 0; l < aas.length; l++) {
						for (int m = 0; m < aas.length; m++) {
							String seq=new String(new char[] {'P', 'K', aas[i], aas[j], aas[k], aas[l], aas[m], 'K'});
							writer.println(seq+","+nce+","+charge);
						}
					}
				}
			}
		}
		writer.close();
	}
	
	public static void main2(String[] args) throws Exception {
		AminoAcidConstants aaConstants=new AminoAcidConstants();
		
		File parent=new File("/Users/searleb/Documents/mak/");
		File fileName=new File("/Users/searleb/Documents/mak/prosit.csv");
		
		float defaultNCE=33.0f;
		byte defaultCharge=3;
		byte[] chargeStates=new byte[] {2,3,4};
		File[] files=new File[] {new File(parent, "metzyme_p2.txt"), new File(parent, "metzyme_p3.txt"), new File(parent, "metzyme_p4.txt")};
		files=new File[] {new File(parent, "ProteOMZ_unique_peptides_exclusive-counts_clean_no-DECOY.csv.txt"),new File(parent, "ProteOMZ_unique_peptides_exclusive-counts_clean_no-DECOY.csv.txt"),new File(parent, "ProteOMZ_unique_peptides_exclusive-counts_clean_no-DECOY.csv.txt")};

		PrintWriter writer=new PrintWriter(fileName);
		writer.println("modified_sequence,collision_energy,precursor_charge");
		int total=0;
		
		for (int i=0; i<files.length; i++) {
			byte charge=chargeStates[i];
			BufferedReader in=new BufferedReader(new FileReader(files[i]));

			String seq;
			while ((seq=in.readLine())!=null) {
				if (seq.trim().length()==0) {
					continue;
				}
				
				if (seq.length()>30||seq.length()<7) continue;
				if (seq.indexOf('B')>=0||seq.indexOf('J')>=0||seq.indexOf('O')>=0||seq.indexOf('U')>=0||seq.indexOf('X')>=0||seq.indexOf('Z')>=0||seq.indexOf('*')>=0) continue;
				
				
				double mz=aaConstants.getChargedMass(seq, charge);
				if (mz>1002.70||mz<396.43) continue;
				
				writer.println(seq+","+convertNCE(defaultNCE, (byte)charge, defaultCharge)+","+(charge));
				total++;
			}
			
			in.close();
		}
		
		writer.close();
		Logger.logLine("Finished writing "+total+" peptides to Prosit CSV!");
	}

	public static void main3(String[] args) throws Exception {
		//File f=new File("/Users/bsearle/Documents/prosit/Pfalciparum/PlasmoDB-43_Pfalciparum3D7_AnnotatedProteins_042419.fasta");
		File fasta=new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta");
		//File f=new File("/Volumes/searle_ssd/malaria/uniprot_human_25apr2019.fasta");
		//File f=new File("/Volumes/searle_ssd/malaria/PlasmoDB-43_Pfalciparum3D7_AnnotatedProteins_042419.fasta");
		//File f=new File("/Users/searleb/Downloads/uniprot-taxonomy_183190.fasta");
		//File fasta=new File("/Users/searleb/Downloads/2019.05_UP000028761_9555_Papio_anubis_canonical_fixed.fasta");
		//fasta=new File("/Users/searleb/Documents/school/xcordia/hela_specific_database/HeLa_Database.txt");
		fasta=new File("/Volumes/searle_ssd/malaria/hela_specific_database/extra_entries.fasta");
		
		int defaultNCE = 33;
		byte defaultCharge = (byte)2;
		byte minCharge=2;
		byte maxCharge=3;
		int maxMissedCleavages=1;
		double minimumMz = 396.43;
		double maximumMz = 1002.70;
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		
		writeCSV(fasta, enzyme, defaultNCE, defaultCharge, minCharge, maxCharge, maxMissedCleavages, new Range(minimumMz, maximumMz), true, false);
	}
	
	public static void main(String[] args) throws Exception {
		File fastaFile=new File("/Users/searleb/Documents/OSU/cptac_grant/cptac_rodland/refseq_human_hg19_refGene_20170329_protein_add_contaminants.fasta");
		File libraryFile=new File("/Users/searleb/Documents/OSU/cptac_grant/cptac_rodland/uniprot_human_25apr2019.z3_nce.dlib");

		LibraryFile library=new LibraryFile();
		library.openFile(libraryFile);
		ArrayList<PeptidePrecursor> entries=library.getAllPeptidePrecursors(new AminoAcidConstants());

		int defaultNCE = 33;
		byte defaultCharge = (byte)3;
		byte minCharge=2;
		byte maxCharge=3;
		int maxMissedCleavages=1;
		double minimumMz = 396.43;
		double maximumMz = 1002.70;
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		HashSet<PeptidePrecursor> peptides=getPeptidesFromFASTA(fastaFile, enzyme, minCharge, maxCharge, maxMissedCleavages, new Range(minimumMz, maximumMz));
		int totalPeptides=peptides.size();
		peptides.removeAll(entries);
		System.out.println("FOUND previous: "+entries.size()+", new: "+totalPeptides+", after removal: "+peptides.size());
		
		String fileName = checkCSVName(null, fastaFile, enzyme, defaultNCE, defaultCharge);
		int total = writePrositFile(fileName, defaultNCE, defaultCharge, true, false, peptides);
		Logger.logLine("Finished writing "+total+" peptides to Prosit CSV!");
	}
	
	public static void writeCSV(File fasta) throws FileNotFoundException {
		int defaultNCE = 33;
		byte defaultCharge = (byte)3;
		byte minCharge=2;
		byte maxCharge=3;
		int maxMissedCleavages=1;
		double minimumMz = 396.43;
		double maximumMz = 1002.70;
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		writeCSV(fasta, enzyme, defaultNCE, defaultCharge, minCharge, maxCharge, maxMissedCleavages, new Range(minimumMz, maximumMz), true, false);
	}

	public static void writeCSV(File fasta, DigestionEnzyme enzyme, int defaultNCE, byte defaultCharge, byte minCharge, byte maxCharge, int maxMissedCleavages, Range mzRange, boolean adjustNCEForDIA, boolean addDecoys) throws FileNotFoundException {
		writeCSV(null, fasta, enzyme, defaultNCE, defaultCharge, minCharge, maxCharge, maxMissedCleavages, mzRange, adjustNCEForDIA, addDecoys);
	}

	public static void writeCSV(String csvFileName, File fasta, DigestionEnzyme enzyme, int defaultNCE, byte defaultCharge, byte minCharge, byte maxCharge, int maxMissedCleavages, Range mzRange, boolean adjustNCEForDIA, boolean addDecoys) throws FileNotFoundException {
		writeCSV(csvFileName, fasta, enzyme, defaultNCE, defaultCharge, minCharge, maxCharge, maxMissedCleavages, mzRange, adjustNCEForDIA, addDecoys, new HashSet<>());
	}
	
	public static void writeCSV(String csvFileName, File fasta, DigestionEnzyme enzyme, int defaultNCE, byte defaultCharge, byte minCharge, byte maxCharge, int maxMissedCleavages, Range mzRange, boolean adjustNCEForDIA, boolean addDecoys, HashSet<PeptidePrecursor> peptidesToIgnore) throws FileNotFoundException {
		String fileName = checkCSVName(csvFileName, fasta, enzyme, defaultNCE, defaultCharge);
		Logger.logLine("Starting to build Prosit CSV: "+fileName);

		HashSet<PeptidePrecursor> allPeptides=getPeptidesFromFASTA(fasta, enzyme, minCharge, maxCharge, maxMissedCleavages, mzRange);
		
		allPeptides.removeAll(peptidesToIgnore);

		int total = writePrositFile(fileName, defaultNCE, defaultCharge, adjustNCEForDIA, addDecoys, allPeptides);
		Logger.logLine("Finished writing "+total+" peptides to Prosit CSV!");
	}

	protected static HashSet<PeptidePrecursor> getPeptidesFromFASTA(File fasta, DigestionEnzyme enzyme, byte minCharge, byte maxCharge, int maxMissedCleavages, Range mzRange) {
		HashSet<PeptidePrecursor> allPeptides=new HashSet<>();
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(fasta, parameters);
		AminoAcidConstants aminoAcidConstants = new AminoAcidConstants();
		
		for (FastaEntryInterface entry : entries) {
			ArrayList<FastaPeptideEntry> peptidesInProtein=enzyme.digestProtein(entry, 7, 30, maxMissedCleavages, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
			for (FastaPeptideEntry pep : peptidesInProtein) {
				for (byte pepCharge = minCharge; pepCharge <=maxCharge; pepCharge++) {
					String seq=pep.getSequence();
					double pepMass=aminoAcidConstants.getMass(seq)+MassConstants.oh2;
					double pepChargedMass=(pepMass+MassConstants.protonMass*pepCharge)/pepCharge;

					if (mzRange.contains(pepChargedMass)) {
						if (seq.indexOf('B')>=0||seq.indexOf('J')>=0||seq.indexOf('O')>=0||seq.indexOf('U')>=0||seq.indexOf('X')>=0||seq.indexOf('Z')>=0||seq.indexOf('*')>=0) {
							continue;
						} else {
							allPeptides.add(new SimplePeptidePrecursor(seq, pepCharge, aminoAcidConstants));
						}
					}
				}
			}
		}
		
		return allPeptides;
	}

	public static void writeCSV(LibraryFile library, int defaultNCE, byte defaultCharge, boolean adjustNCEForDIA, boolean addDecoys) throws FileNotFoundException {
		writeCSV(null, library, defaultNCE, defaultCharge, adjustNCEForDIA, addDecoys);
	}
	public static void writeCSV(String csvFileName, LibraryFile library, int defaultNCE, byte defaultCharge, boolean adjustNCEForDIA, boolean addDecoys) throws FileNotFoundException {
		String fileName = checkCSVName(csvFileName, library.getFile(), null, defaultNCE, defaultCharge);
		Logger.logLine("Starting to build Prosit CSV: "+fileName);
		AminoAcidConstants aminoAcidConstants = new AminoAcidConstants();

		HashSet<PeptidePrecursor> allPeptides=new HashSet<>();
		try {
			ArrayList<LibraryEntry> allEntries = library.getAllEntries(false, aminoAcidConstants);
			allPeptides.addAll(allEntries);
		} catch (IOException e) {
			throw new EncyclopediaException("Error parsing library", e);
		} catch (SQLException e) {
			throw new EncyclopediaException("Error parsing library", e);
		} catch (DataFormatException e) {
			throw new EncyclopediaException("Error parsing library", e);
		}
		
		int total = writePrositFile(fileName, defaultNCE, defaultCharge, adjustNCEForDIA, addDecoys, allPeptides);
		Logger.logLine("Finished writing "+total+" peptides to Prosit CSV!");
	}

	/**
	 * 
	 * @param csvFileName CAN BE NULL (if null, then create a name)
	 * @param fasta
	 * @param enzyme CAN BE NULL (if null, then use blank text in the created name)
	 * @param defaultNCE
	 * @param defaultCharge
	 * @return
	 */
	static String checkCSVName(String csvFileName, File fasta, DigestionEnzyme enzyme, int defaultNCE, byte defaultCharge) {
		String fileName;
		if (null==csvFileName||StringUtils.isBlank(csvFileName)) {
			String enzymeText=enzyme==null?"":("."+enzyme.getPercolatorName());
			fileName = fasta.getAbsolutePath() + enzymeText + ".z" + defaultCharge + "_nce" + defaultNCE + ".csv";
		} else {
			fileName = csvFileName;
		}
		return fileName;
	}
	
	public static int writePrositFile(String fileName, int defaultNCE, byte defaultCharge, boolean adjustNCEForDIA, boolean addDecoys, HashSet<PeptidePrecursor> allPeptides) throws FileNotFoundException {
		AminoAcidConstants constants = new AminoAcidConstants(); // does not support PTMs
		PrintWriter writer=new PrintWriter(fileName);
		int total=0;
		HashSet<PeptidePrecursor> writtenPeptides=new HashSet<>();
		
		writer.println("modified_sequence,collision_energy,precursor_charge");
		for (PeptidePrecursor peptidePrecursor : allPeptides) {
			
			// Prosit only supports unmodified peptides:
			String seq = peptidePrecursor.getPeptideSeq();
			byte precursorCharge = peptidePrecursor.getPrecursorCharge();
			PeptidePrecursor unmodified=new SimplePeptidePrecursor(seq, precursorCharge, constants);
			if (writtenPeptides.contains(unmodified)) {
				continue;
			}

			// Prosit doesn't support charge >6
			if (precursorCharge<1&&precursorCharge>6) {
				continue;
			}
			
			// remove peptides that don't match PROSIT limitations:
			if (seq.indexOf('B')>=0||seq.indexOf('J')>=0||seq.indexOf('O')>=0||seq.indexOf('U')>=0||seq.indexOf('X')>=0||seq.indexOf('Z')>=0||seq.indexOf('*')>=0) {
				continue;
			}
			if (seq.length()<7||seq.length()>30) {
				continue;
			}
			
			writtenPeptides.add(unmodified);
			
			float nce = adjustNCEForDIA?convertNCE(defaultNCE, precursorCharge, defaultCharge):defaultNCE;
			writer.println(seq+","+nce+","+precursorCharge);
			if (addDecoys) {
				String reverse=PeptideUtils.reverse(seq, constants);
				writer.println(reverse+","+nce+","+precursorCharge);
			}
			total++;
		}
		writer.close();
		return total;
	}
	
	// http://proteomicsnews.blogspot.com/2014/06/normalized-collision-energy-calculation.html
	public static float convertNCE(float nce, byte charge, byte defaultCharge) {
		return nce*getChargeFactor(defaultCharge)/getChargeFactor(charge);
	}

	private static float getChargeFactor(byte charge) {
		switch (charge) {
			case 1: return 1.0f;
			case 2: return 0.9f;
			case 3: return 0.85f;
			case 4: return 0.8f;
			case 5: return 0.75f;
			default: return 0.75f;
		}
	}
}
