package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Correlation;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class LibraryComparisonTest {
	public static void main(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();

//		File[] libraryFilesDDA=new File[] {
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce15.dlib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce18.dlib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce21.dlib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce24.dlib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce27.dlib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce30.dlib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce33.dlib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce36.dlib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce39.dlib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce42.dlib"),
//				new File("/Volumes/searle_ssd/malaria/DDA_yeast_with_iRTs.dlib")
//		};
//		File[] libraryFiles=new File[] {
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce15_clib.elib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce18_clib.elib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce21_clib.elib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce24_clib.elib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce27_clib.elib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce30_clib.elib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce33_clib.elib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce36_clib.elib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce39_clib.elib"),
//				new File("/Volumes/searle_ssd/malaria/uniprot_yeast_25jan2019.fasta.z2_nce42_clib.elib"),
//				new File("/Volumes/searle_ssd/malaria/DDA_yeast_with_iRTs.dlib")
//		};
		
		File[] libraryFilesDDA=new File[] {
				new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce15.dlib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce18.dlib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce21.dlib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce24.dlib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce27.dlib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce30.dlib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce33.dlib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce36.dlib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce39.dlib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce42.dlib"),
				//new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/versus_high_phRP/yeast_scx_dda_tpp.dlib")
		};
		File[] libraryFilesDIA=new File[] {
				new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce15_clib.elib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce18_clib.elib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce21_clib.elib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce24_clib.elib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce27_clib.elib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce30_clib.elib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce33_clib.elib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce36_clib.elib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce39_clib.elib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce42_clib.elib"),
				//new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/versus_high_phRP/yeast_scx_dda_tpp.dlib")
		};
		
		File[] libraryFiles=libraryFilesDIA;
		
		File diaFile=new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/raw_files/20190206_LUM1_CPBA_EASY04_060_30_SA_90mingrad_80B_DIA_400_1000_8mzol_15k_20IIT_4e5agc_1633-01_01.mzML.elib");
		
		LibraryFile[] libraries=new LibraryFile[libraryFiles.length];
		for (int i=0; i<libraries.length; i++) {
			libraries[i]=new LibraryFile();
			libraries[i].openFile(libraryFiles[i]);
		}
		
		LibraryFile file=new LibraryFile();
		file.openFile(diaFile);
		ArrayList<LibraryEntry> entries=file.getAllEntries(false, parameters.getAAConstants());
		
		System.out.println("Processing "+entries.size());
		System.out.print("count,peptide,charge,protein");
		for (int i=0; i<libraries.length; i++) {
			System.out.print(","+libraryFiles[i].getName());
		}
		System.out.println();

		int count=0;
		for (LibraryEntry entry : entries) {
			if (entry.getMassArray().length<6||General.sum(entry.getIntensityArray())<=0.0f) continue;
			
			count++;
			StringBuilder sb=new StringBuilder(count+","+entry.getPeptideModSeq()+","+entry.getPrecursorCharge()+","+General.toString(entry.getAccessions().toArray(), ";"));

			boolean skip=Math.random()>0.05f;
			for (int i=0; i<libraries.length; i++) {
				if (!skip) {
					ArrayList<LibraryEntry> candidates=libraries[i].getEntries(entry.getPeptideModSeq(), entry.getPrecursorCharge(), false);
					if (candidates.size()>0) {
						AnnotatedLibraryEntry dda=AnnotatedLibraryEntry.getAnnotationsOnly(candidates.get(0), parameters);
						float correlation=(float)Correlation.getPearsons(entry, dda, parameters.getFragmentTolerance());
						sb.append(","+correlation);
					} else {
						sb.append(",-2");
						//skip=true;
					}
				}
			}
			if (!skip) {
				System.out.println(sb.toString());
			}
		}
		file.close();
	}
}
