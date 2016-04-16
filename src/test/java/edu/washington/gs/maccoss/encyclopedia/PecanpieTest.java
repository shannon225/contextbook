package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;

public class PecanpieTest {
	public static void main(String[] args) {
		// EXAMPLE
		File diaFile=new File("/Users/searleb/Documents/school/projects/pecandata/DIA_1xGFP_20x20mz_500to900_rep1.mzML");
		File fastaFile=new File("/Users/searleb/Documents/school/projects/pecandata/UP000005640_9606.fasta");
		File featureFile=new File("/Users/searleb/Documents/school/projects/pecandata/encyc_report.feature.txt");
		File outputFile=new File("/Users/searleb/Documents/school/projects/pecandata/encyc_report.percolator.txt");
		PecanSearchParameters parameters=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.YONLY, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"));
		PecanScoringFactory factory=new PecanOneScoringFactory(parameters, featureFile);
		
		ArrayList<FastaPeptideEntry> targets=new ArrayList<FastaPeptideEntry>();
		targets.add(new FastaPeptideEntry("FILE", ">Protein", "IGHTVEREDTPAIR"));
		targets=null;
		
		try {
			Pecanpie.runPie(new EmptyProgressIndicator(), Optional.ofNullable(targets), diaFile, fastaFile, featureFile, outputFile, factory);
		} catch (Exception e) {
			System.err.println("Encountered Fatal Error!");
			e.printStackTrace();
		}
	}
}
