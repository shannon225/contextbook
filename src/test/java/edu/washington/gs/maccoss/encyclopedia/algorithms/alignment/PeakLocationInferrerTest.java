package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;

public class PeakLocationInferrerTest {
	public static void main(String[] args) throws Exception {
		//LibraryFile.OPEN_IN_PLACE=true;
		
		HashMap<String, String> map=SearchParameterParser.getDefaultParameters();
		SearchParameters parameters=SearchParameterParser.parseParameters(map);
		
		//QuantitativeSearchJobData job1=getData(parameters, "/Users/searleb/Documents/projects/encyclopedia/mzml/zero_hela/121115_bcs_hela_24mz_400_1000_0D_1.dia");
		//QuantitativeSearchJobData job2=getData(parameters, "/Users/searleb/Documents/projects/encyclopedia/mzml/zero_hela/121115_bcs_hela_24mz_400_1000_0D_2.dia");
		QuantitativeSearchJobData job1=getData(parameters, "/Users/searleb/Documents/school/projects/may_asms/hela/on_column/timecourse/23aug2017_hela_serum_timecourse_wide_1a.dia");
		QuantitativeSearchJobData job2=getData(parameters, "/Users/searleb/Documents/school/projects/may_asms/hela/on_column/timecourse/23aug2017_hela_serum_timecourse_wide_1b.dia");
		ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
		jobs.add(job1);
		jobs.add(job2);
		
		PeakLocationInferrerInterface inferrer=AlternatePeakLocationInferrer.getAlignmentData(new EmptyProgressIndicator(), jobs, getPeptides(jobs), parameters);
		System.out.println("both j1: "+inferrer.getPreciseRTInSec(job1, "GADGMILGGPQSDSDTDAQR", 3085f)+"\t GADGMILGGPQSDSDTDAQR"); // both
		System.out.println("both j2: "+inferrer.getPreciseRTInSec(job2, "GADGMILGGPQSDSDTDAQR", 3085f)+"\t GADGMILGGPQSDSDTDAQR"); // both
		System.out.println("a    j1: "+inferrer.getPreciseRTInSec(job1, "GPPAPTTQAQPDLIKPLPLHK", 3517f)+"\t GPPAPTTQAQPDLIKPLPLHK"); // only in a
		System.out.println("a    j2: "+inferrer.getPreciseRTInSec(job2, "GPPAPTTQAQPDLIKPLPLHK", 3517f)+"\t GPPAPTTQAQPDLIKPLPLHK"); // only in a
		System.out.println("b    j1: "+inferrer.getPreciseRTInSec(job1, "NSSYVHGGVDASGKPQEAVYGQNDIHHK", 2198f)+"\t NSSYVHGGVDASGKPQEAVYGQNDIHHK"); // only in b
		System.out.println("b    j2: "+inferrer.getPreciseRTInSec(job2, "NSSYVHGGVDASGKPQEAVYGQNDIHHK", 2198f)+"\t NSSYVHGGVDASGKPQEAVYGQNDIHHK"); // only in b	
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
	
	private static ArrayList<PercolatorPeptide> getPeptides(ArrayList<SearchJobData> jobs) {
		ArrayList<PercolatorPeptide> peptides=new ArrayList<>();
		
		for (SearchJobData job : jobs) {
			ArrayList<PercolatorPeptide> local=PercolatorReader.getPassingPeptidesFromTSV(job.getOutputFile(), job.getParameters().getEffectivePercolatorThreshold(), false);
			for (PercolatorPeptide pep : local) {
				if ("KGSITSVQAIYVPADDLTDPAPATTFAHLDATTVLSR".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("NSSYVHGGVDASGKPQEAVYGQNDIHHK".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("YYIQNGIQSFMQNYSSIDVLLHQSR".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("HAVSDPSILDSLDLNEDEREVLINNINRR".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("ERVEAVNMAEGIIHDTETK".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("FLNEHPGGEEVLLEQAGVDASESFEDVGHSSDAR".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("TDQVIQSLIALVNDPQPEHPLR".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("TFSHELSDFGLESTAGEIPVVAIR".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("KLEDQLQGGQLEEVILQAEHELNLAR".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("TLIENGEK".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("GGVDVTLPR".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("TVGVEPAADGK".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("IKGDVDVSVPEVEGK".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("VTFEELR".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("GADGMILGGPQSDSDTDAQR".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("HRGQAAQPEPSTGFTATPPAPDSPQEPLVLR".equals(pep.getPeptideModSeq())) peptides.add(pep);
				else if ("GPPAPTTQAQPDLIKPLPLHK".equals(pep.getPeptideModSeq())) peptides.add(pep);
				
				
			}	
		}
		Collections.sort(peptides, new Comparator<PercolatorPeptide>() {
			@Override
			public int compare(PercolatorPeptide o1, PercolatorPeptide o2) {
				if (o1==null&&o2==null) return 0;
				if (o1==null) return -1;
				if (o2==null) return 1;
				
				return o1.getPeptideModSeq().compareTo(o2.getPeptideModSeq());
			}
		});
		for (PercolatorPeptide pep : peptides) {
			System.out.println("found "+pep.getPsmID());
		}
		
		return peptides;
	}
}
