package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import junit.framework.TestCase;

public class PeakLocationInferrerTest {
	public static void main(String[] args) throws Exception {
		LibraryFile.OPEN_IN_PLACE=true;
		
		HashMap<String, String> map=SearchParameterParser.getDefaultParameters();
		SearchParameters parameters=SearchParameterParser.parseParameters(map);
		
		QuantitativeSearchJobData job1=getData(parameters, "/Users/searleb/Documents/projects/encyclopedia/mzml/zero_hela/121115_bcs_hela_24mz_400_1000_0D_1.dia");
		QuantitativeSearchJobData job2=getData(parameters, "/Users/searleb/Documents/projects/encyclopedia/mzml/zero_hela/121115_bcs_hela_24mz_400_1000_0D_2.dia");
		
		ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
		jobs.add(job1);
		jobs.add(job2);
		
		PeakLocationInferrer inferrer=PeakLocationInferrer.getAlignmentData(new EmptyProgressIndicator(), jobs, getPeptides(), parameters);
		System.out.println("j1: "+inferrer.getPreciseRTInSec(job1, "NSSYVHGGVDASGKPQEAVYGQNDIHHK", 2528f));
		System.out.println("j2: "+inferrer.getPreciseRTInSec(job2, "NSSYVHGGVDASGKPQEAVYGQNDIHHK", 2528f));
	}

	private static QuantitativeSearchJobData getData(SearchParameters parameters, String dia) {
		File diaFile=new File(dia);
		File outputFile=new File(diaFile.getAbsolutePath()+".encyclopedia.txt");
		File decoyFile=new File(diaFile.getAbsolutePath()+".encyclopedia.decoy.txt");
		File featureFile=new File(diaFile.getAbsolutePath()+".features.txt");
		final File libraryFile=new File(diaFile.getAbsolutePath()+".elib");
		QuantitativeSearchJobData job=new QuantitativeSearchJobData(diaFile, featureFile, outputFile, decoyFile, parameters, "1") {
			@Override
			public String getSearchType() {
				return "Dummy";
			}
			
			@Override
			public File getResultLibrary() {
				return libraryFile;
			}
		};
		return job;
	}
	
	private static ArrayList<PercolatorPeptide> getPeptides() {
		ArrayList<PercolatorPeptide> peptides=new ArrayList<>();
		peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_1.dia:4922.879:KGSITSVQAIYVPADDLTDPAPATTFAHLDATTVLSR+4"));
		//peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_2.dia:4936.5684:KGSITSVQAIYVPADDLTDPAPATTFAHLDATTVLSR+4"));
		
		peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_1.dia:2528.0693:NSSYVHGGVDASGKPQEAVYGQNDIHHK+3"));
		//peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_2.dia:2557.912:NSSYVHGGVDASGKPQEAVYGQNDIHHK+3"));
		
		//peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_1.dia:5767.299:YYIQNGIQSFMQNYSSIDVLLHQSR+3"));
		peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_2.dia:5776.2207:YYIQNGIQSFMQNYSSIDVLLHQSR+3"));
		
		//peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_1.dia:4737.8193:HAVSDPSILDSLDLNEDEREVLINNINRR+4"));
		peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_2.dia:4750.6094:HAVSDPSILDSLDLNEDEREVLINNINRR+4"));
		
		peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_1.dia:4225.2847:ERVEAVNMAEGIIHDTETK+3"));
		//peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_2.dia:4232.4346:ERVEAVNMAEGIIHDTETK+3"));
		
		peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_1.dia:4915.851:FLNEHPGGEEVLLEQAGVDASESFEDVGHSSDAR+4"));
		//peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_2.dia:4922.305:FLNEHPGGEEVLLEQAGVDASESFEDVGHSSDAR+4"));
		
		//peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_1.dia:5556.0767:TDQVIQSLIALVNDPQPEHPLR+3"));
		peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_2.dia:5566.2466:TDQVIQSLIALVNDPQPEHPLR+3"));
		
		//peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_1.dia:5013.528:TFSHELSDFGLESTAGEIPVVAIR+3"));
		peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_2.dia:5027.694:TFSHELSDFGLESTAGEIPVVAIR+3"));
		
		//peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_1.dia:5782.0054:KLEDQLQGGQLEEVILQAEHELNLAR+3"));
		peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_2.dia:5796.0913:KLEDQLQGGQLEEVILQAEHELNLAR+3"));
		
		peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_2.dia:2256.097:TLIENGEK+2"));
		peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_2.dia:3308.137:GGVDVTLPR+2"));
		
		peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_1.dia:2157.279:TVGVEPAADGK+2"));
		peptides.add(new DummyPercolatorPeptide("121115_bcs_hela_24mz_400_1000_0D_1.dia:3238.2861:IKGDVDVSVPEVEGK+2"));
		
		return peptides;
	}
	
	private static class DummyPercolatorPeptide extends PercolatorPeptide {
		public DummyPercolatorPeptide(String psmID) {
			super(psmID, "p1", 0.001f, 0.01f);
		}
		
	}
}
