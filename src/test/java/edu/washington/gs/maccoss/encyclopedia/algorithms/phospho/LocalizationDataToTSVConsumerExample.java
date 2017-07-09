package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;

public class LocalizationDataToTSVConsumerExample {
	public static void main(String[] args) {
		File f=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/22jun2016_mcf7_phospho_1a.dia.thesaurus.txt.localizations.txt");
		File p=new File("/Users/searleb/Documents/school/localization_manuscript/mcf7/22jun2016_mcf7_phospho_1a.dia.thesaurus.txt");
		ArrayList<PercolatorPeptide> passingPeptides=PercolatorReader.getPassingPeptidesFromTSV(p, 0.05f);
		System.out.println("Found "+passingPeptides.size()+" total peptides...");

		HashMap<String, String> defaults=SearchParameterParser.getDefaultParameters();
		defaults.put("-localizationModification", "Phosphorylation");
		defaults.put("-ptol", "16.67");
		defaults.put("-ftol", "16.67");
		defaults.put("-lftol", "16.67");
		//defaults.put("-frag", "yonly");
		defaults.put("-scoringBreadthType", "uncal20");
		SearchParameters parameters=SearchParameterParser.parseParameters(defaults);
		
		HashMap<String, ModificationLocalizationData> localizationData=LocalizationDataToTSVConsumer.readLocalizationFile(f, passingPeptides, parameters);
		System.out.println("Found "+localizationData.size()+" localized peptides...");
	}

}
