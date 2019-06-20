package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class PrositCSVWriter {

	public static void main(String[] args) throws Exception {
		//File f=new File("/Users/bsearle/Documents/prosit/Pfalciparum/PlasmoDB-43_Pfalciparum3D7_AnnotatedProteins_042419.fasta");
		//File f=new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta");
		//File f=new File("/Volumes/searle_ssd/malaria/uniprot_human_25apr2019.fasta");
		//File f=new File("/Volumes/searle_ssd/malaria/PlasmoDB-43_Pfalciparum3D7_AnnotatedProteins_042419.fasta");
		//File f=new File("/Users/searleb/Downloads/uniprot-taxonomy_183190.fasta");
		File fasta=new File("/Users/searleb/Downloads/2019.05_UP000028761_9555_Papio_anubis_canonical_fixed.fasta");
		
		int defaultNCE = 33;
		byte defaultCharge = (byte)3;
		int minCharge=2;
		int maxCharge=3;
		int maxMissedCleavages=1;
		
		writeCSV(fasta, defaultNCE, defaultCharge, minCharge, maxCharge, maxMissedCleavages);
	}
	
	public static void writeCSV(File fasta) throws FileNotFoundException {
		int defaultNCE = 33;
		byte defaultCharge = (byte)3;
		int minCharge=2;
		int maxCharge=3;
		int maxMissedCleavages=1;
		writeCSV(fasta, defaultNCE, defaultCharge, minCharge, maxCharge, maxMissedCleavages);
	}

	public static void writeCSV(File fasta, int defaultNCE, byte defaultCharge, int minCharge, int maxCharge, int maxMissedCleavages) throws FileNotFoundException {
		int[] chargeStates = new int[maxCharge-minCharge+1];
		for (int i = 0; i < chargeStates.length; i++) {
			chargeStates[i]=i+minCharge;
			Logger.logLine("Considering charge +"+chargeStates[i]+"H");
		}
		
		String fileName = fasta.getAbsolutePath()+".z"+defaultCharge+"_nce"+defaultNCE+".csv";
		PrintWriter writer=new PrintWriter(fileName);
		Logger.logLine("Starting to build Prosit CSV: "+fileName);

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(fasta, parameters);
		AminoAcidConstants constants=new AminoAcidConstants();
		
		HashSet<String>[] allPeptides=new HashSet[chargeStates.length];
		for (int i=0; i<allPeptides.length; i++) {
			allPeptides[i]=new HashSet<>();
		}
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		for (FastaEntryInterface entry : entries) {
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 7, 30, maxMissedCleavages, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
			for (FastaPeptideEntry pep : peptides) {
				for (int pepCharge : chargeStates) {
					String seq=pep.getSequence();
					double pepMass=constants.getMass(seq)+MassConstants.oh2;
					double pepChargedMass=(pepMass+MassConstants.protonMass*pepCharge)/pepCharge;

					if (pepChargedMass>(396.43)&&pepChargedMass<(1002.70)) {
						if (seq.indexOf('B')>=0||seq.indexOf('J')>=0||seq.indexOf('O')>=0||seq.indexOf('U')>=0||seq.indexOf('X')>=0||seq.indexOf('Z')>=0||seq.indexOf('*')>=0) {
							continue;
						} else {
							allPeptides[pepCharge-chargeStates[0]].add(seq);
						}
					}
				}
			}
		}
		
		int total=0;
		writer.println("modified_sequence,collision_energy,precursor_charge");
		for (int i=0; i<allPeptides.length; i++) {
			total+=allPeptides[i].size();
			int charge=i+minCharge;
			for (String string : allPeptides[i]) {
				writer.println(string+","+convertNCE(defaultNCE, (byte)charge, defaultCharge)+","+(charge));
			}
		}
		writer.close();
		Logger.logLine("Finished writing "+total+" peptides to Prosit CSV!");
	}
	
	// http://proteomicsnews.blogspot.com/2014/06/normalized-collision-energy-calculation.html
	private static float convertNCE(float nce, byte charge, byte defaultCharge) {
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
