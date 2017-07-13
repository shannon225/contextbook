package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.lang3.StringUtils;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.map.hash.TObjectIntHashMap;

public class LocalizationDataToTSVConsumerExample {
	public static void main(String[] args) {
		//File f=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/22jun2016_mcf7_phospho_1a.dia.thesaurus.txt.localizations.txt");
		//File p=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/22jun2016_mcf7_phospho_1a.dia.thesaurus.txt");

		//File f=new File("/Users/searleb/Documents/school/localization_manuscript/phospho_repeats/01_error/20170430_HeLa_phosp_DIA_B_01_170506220515.mzML.thesaurus.txt.localizations.txt");
		//File p=new File("/Users/searleb/Documents/school/localization_manuscript/phospho_repeats/01_error/20170430_HeLa_phosp_DIA_B_01_170506220515.mzML.thesaurus.txt");

		File[] fs=new File("/Users/searleb/Documents/school/localization_manuscript/phospho_repeats/final_data/hela_repeats").listFiles();
		for (File p : fs) {
			if (p.getName().endsWith(".thesaurus.txt")) {

				File f=new File(p.getAbsolutePath()+".localizations.txt");
				
				ArrayList<PercolatorPeptide> passingPeptides=PercolatorReader.getPassingPeptidesFromTSV(p, 0.05f, false);
				HashSet<String> peptides=new HashSet<>();
				for (PercolatorPeptide peptide : passingPeptides) {
					int numMods=StringUtils.countMatches(peptide.getPeptideModSeq(), "[");
					peptides.add(numMods+PeptideUtils.getPeptideSeq(peptide.getPeptideModSeq()));
				}
				
				System.out.println(p.getName()+" Found "+passingPeptides.size()+" total peptides ("+peptides.size()+" unique)...");

				HashMap<String, String> defaults=SearchParameterParser.getDefaultParameters();
				defaults.put("-localizationModification", "Phosphorylation");
				defaults.put("-ptol", "16.67");
				defaults.put("-ftol", "16.67");
				defaults.put("-lftol", "16.67");
				//defaults.put("-frag", "yonly");
				defaults.put("-scoringBreadthType", "uncal20");
				SearchParameters parameters=SearchParameterParser.parseParameters(defaults);
				
				HashMap<String, ModificationLocalizationData> localizationData=LocalizationDataToTSVConsumer.readLocalizationFile(f, passingPeptides, parameters);
				System.out.println("Found "+localizationData.size()+" peptides in localization file...");
				
				peptides=new HashSet<>();
				int count=0;
				for (ModificationLocalizationData data : localizationData.values()) {
					if (data.isSiteSpecific()) {
						count++;
						peptides.add(data.getNumberOfMods()+PeptideUtils.getPeptideSeq(data.getLocalizationPeptideModSeq().getPeptideModSeq()));
					}
				}
				System.out.println("Found "+count+" site localized peptides ("+peptides.size()+" unique)...");
				
				TObjectIntHashMap<String> multipleForms=new TObjectIntHashMap<>();
				for (ModificationLocalizationData data : localizationData.values()) {
					if (data.isSiteSpecific()&&data.getNumberOfMods()==1&&data.getLocalizationPeptideModSeq().getNumModifiableSites()>1) {
						multipleForms.adjustOrPutValue(data.getNumberOfMods()+PeptideUtils.getPeptideSeq(data.getLocalizationPeptideModSeq().getPeptideModSeq()), 1, 1);
					}
				}
				
				int[] values=multipleForms.values();
				int[] counts=new int[10];
				for (int i=0; i<values.length; i++) {
					counts[values[i]]++;
				}
				System.out.println("Singly Modified Form Histogram:");
				for (int i=0; i<counts.length; i++) {
					System.out.println(i+"\t"+counts[i]);
				}
			}
		}
		
	}

}
