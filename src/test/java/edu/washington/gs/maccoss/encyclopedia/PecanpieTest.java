package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.util.ArrayList;

import com.google.common.base.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;

public class PecanpieTest {
	public static void main(String[] args) {
		// EXAMPLE
		File diaFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/20150708_Ecoli_0911_25x4mzDIA_500_600.dia");
		//File fastaFile=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/ecoli-190209-contam_correctNL.fasta");
		File fastaFile=new File("/Users/searleb/Documents/projects/pecan/v0.9.7/ecoli_20150911_uniprot_sp_digested_Mass600to4000.fasta"); //FIXME
		File featureFile=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/encyc_report.feature.txt");
		File outputFile=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/encyc_report.percolator.txt");
		SearchParameters parameters=new SearchParameters(new AminoAcidConstants(), FragmentationType.YONLY, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"));
		PecanScoringFactory factory=new PecanOneScoringFactory(parameters, featureFile);
		
		ArrayList<FastaEntry> targets=new ArrayList<FastaEntry>();
		targets.add(new FastaEntry("FILE", ">Protein", "IGHTVEREDTPAIR"));
		targets=null;
		
		try {
			Pecanpie.runPie(Optional.fromNullable(targets), diaFile, fastaFile, featureFile, outputFile, factory);
		} catch (Exception e) {
			System.err.println("Encountered Fatal Error!");
			e.printStackTrace();
		}
	}
}
