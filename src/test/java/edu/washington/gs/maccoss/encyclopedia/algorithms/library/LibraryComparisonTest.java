package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Ion;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Correlation;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class LibraryComparisonTest {
	public static void main2(String[] args) throws Exception {
		File[] libraryFilesDIA=new File[] {
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.chymotrypsin.z3_nce27.dlib"),
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.chymotrypsin.z3_nce30.dlib"),
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.chymotrypsin.z3_nce33.dlib"),
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.chymotrypsin.z3_nce36.dlib"),
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.chymotrypsin.z3_nce39.dlib"),
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.glu-c.z3_nce27.dlib"),
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.glu-c.z3_nce30.dlib"),
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.glu-c.z3_nce33.dlib"),
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.glu-c.z3_nce36.dlib"),
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.glu-c.z3_nce39.dlib"),
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.trypsin.z3_nce27.dlib"),
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.trypsin.z3_nce30.dlib"),
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.trypsin.z3_nce33.dlib"),
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.trypsin.z3_nce36.dlib"),
				new File("/Users/searleb/Documents/swaney/CMV_CE_prosit_predictions/uniprot-proteome_up000000938.fasta.trypsin.z3_nce39.dlib")};

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		float[][] matrix=new float[libraryFilesDIA.length][];
		for (int f=0; f<libraryFilesDIA.length; f++) {
			int[] histogram=new int[21];
			LibraryFile lib=new LibraryFile();
			lib.openFile(libraryFilesDIA[f]);
			TFloatArrayList counts=new TFloatArrayList();
			for (LibraryEntry entry : lib.getAllEntries(false, parameters.getAAConstants())) {
				AnnotatedLibraryEntry annotated=new AnnotatedLibraryEntry(entry, parameters);
				Ion[] ions=annotated.getIonAnnotations();
				
				float[] intensities=entry.getIntensityArray();
				float target=0.25f*General.max(intensities);
				
				TFloatArrayList altIntensities=new TFloatArrayList();
				for (int i=0; i<ions.length; i++) {
					if (ions[i]!=null) {
						altIntensities.add(intensities[i]);
					}
				}
				intensities=altIntensities.toArray();
				
				
				int countAboveTarget=0;
				for (int i=0; i<intensities.length; i++) {
					if (intensities[i]>target) {
						countAboveTarget++;
					}
				}
				if (histogram.length<=countAboveTarget) {
					histogram[histogram.length-1]++;
				} else {
					histogram[countAboveTarget]++;
				}
				counts.add(countAboveTarget);
			}
			float histogramSum=General.sum(histogram);
			matrix[f]=General.divide(General.toFloatArray(histogram), histogramSum);
			
			float[] array=counts.toArray();
			System.out.println(libraryFilesDIA[f].getName()+" \t"+General.mean(array));
		}

		System.out.print("n");
		for (int f=0; f<libraryFilesDIA.length; f++) {
			System.out.print("\t"+libraryFilesDIA[f].getName());
		}
		System.out.println();
		for (int i=0; i<matrix[0].length; i++) {
			System.out.print(i);
			for (int f=0; f<libraryFilesDIA.length; f++) {
				System.out.print("\t"+matrix[f][i]);
			}
			System.out.println();
		}
	}
	
	public static void main3(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		System.out.println("reading prosit...");
		File prosit=new File("/Volumes/searle_ssd/malaria/hela_specific_database/HeLa_Database.txt.z2_nce33.dlib");
		LibraryFile pfile=new LibraryFile();
		pfile.openFile(prosit);
		ArrayList<LibraryEntry> entries=pfile.getAllEntries(false, parameters.getAAConstants());
		TObjectFloatHashMap<String> rtMap=new TObjectFloatHashMap<>();
		int count=0;
		for (LibraryEntry entry : entries) {
			count++;
			if (count%100000==0) System.out.println(count+"...");
			rtMap.put(entry.getPeptideModSeq(), entry.getRetentionTime());
		}
		pfile.close();
		
		System.out.println("reading massive...");
		File kb=new File("/Volumes/searle_ssd/malaria/MassIVE-KB/LIBRARY_CREATION_AUGMENT_LIBRARY_TEST-82c0124b-download_filtered_sptxt_library-main.dlib");
		LibraryFile kbfile=new LibraryFile();
		kbfile.openFile(kb);
		
		entries=kbfile.getAllEntries(false, parameters.getAAConstants());

		count=0;
		ArrayList<LibraryEntry> calibratedEntries=new ArrayList<>();
		for (LibraryEntry entry : entries) {
			count++;
			if (count%100000==0) System.out.println(count+"...");
			if (rtMap.contains(entry.getPeptideModSeq())) {
				calibratedEntries.add(entry.updateRetentionTime(rtMap.get(entry.getPeptideModSeq())));
			}
		}
		kbfile.close();
		
		System.out.println("writing crosslib...");
		LibraryFile savefile=new LibraryFile();
		savefile.openFile();
		savefile.dropIndices();
		savefile.addEntries(calibratedEntries);
		savefile.createIndices();
		
		savefile.saveAsFile(new File(kb.getParentFile(), "massive_kb_with_prosit_rts.dlib"));
		System.out.println("done!");
	}
	public static void main4(String[] args) throws Exception {
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
				//new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/versus_high_phRP/yeast_hpHRP_dda_tpp.dlib")
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
				new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/versus_high_phRP/yeast_hpHRP_dda_tpp.dlib")
		};
		File[] libraryFileFromUW=new File[] {
				//new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce33.dlib"),
				new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/versus_high_phRP/yeast_hpHRP_dda_tpp.dlib"),
				//new File("/Users/searleb/Downloads/23aug2017_yeast_timecourse_clib.elib")
		};
		File[] libraryFileFromPS=new File[] {
				new File("/Users/searleb/Downloads/myPrositLib_2019.dlib"),
				new File("/Users/searleb/Downloads/myPrositLib_2020.dlib"),
		};
		
		File[] libraryFiles=libraryFileFromPS;
		
		File comparisonFile=new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/raw_files/20190206_LUM1_CPBA_EASY04_060_30_SA_90mingrad_80B_DIA_400_1000_8mzol_15k_20IIT_4e5agc_1633-01_01.mzML.elib");
		comparisonFile=new File("/Volumes/searle_ssd/malaria/novo_yeast/DIA_analysis/clibs_vs_predicted/uniprot_yeast_25jan2019.fasta.z2_nce33_clib.elib");
		comparisonFile=new File("/Users/searleb/Downloads/filtered_mustela_final_merged_with_coronavirus_results.dlib");
		
		LibraryFile[] libraries=new LibraryFile[libraryFiles.length];
		for (int i=0; i<libraries.length; i++) {
			libraries[i]=new LibraryFile();
			libraries[i].openFile(libraryFiles[i]);
		}
		
		LibraryFile file=new LibraryFile();
		file.openFile(comparisonFile);
		ArrayList<LibraryEntry> entries=file.getAllEntries(false, parameters.getAAConstants());
		
		System.out.println("Processing "+entries.size());
		System.out.print("count,peptide,charge,protein");
		for (int i=0; i<libraries.length; i++) {
			System.out.print(","+libraryFiles[i].getName());
		}
		System.out.println();

		int count=0;
		for (LibraryEntry entry : entries) {
			entry=AnnotatedLibraryEntry.getAnnotationsOnly(entry, parameters);
			
			if (entry.getMassArray().length<6||General.sum(entry.getIntensityArray())<=0.0f) continue;
			
			count++;
			StringBuilder sb=new StringBuilder(count+","+entry.getPeptideModSeq()+","+entry.getPrecursorCharge()+","+General.toString(entry.getAccessions().toArray(), ";"));

			boolean skip=false;//Math.random()>0.05f;
			for (int i=0; i<libraries.length; i++) {
				if (!skip) {
					ArrayList<LibraryEntry> candidates=libraries[i].getEntries(entry.getPeptideModSeq(), entry.getPrecursorCharge(), false);
					if (candidates.size()>0) {
						AnnotatedLibraryEntry dda=AnnotatedLibraryEntry.getAnnotationsOnly(candidates.get(0), parameters);
						float correlation=(float)Correlation.getPearsons(entry, dda, parameters.getFragmentTolerance());
						sb.append(","+correlation);
					} else {
						sb.append(",-2");
						skip=true;
					}
				}
			}
			if (!skip) {
				System.out.println(sb.toString());
			}
		}
		file.close();
	}

	public static void main(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		File larger=new File("/Users/searleb/Documents/damien/dda_library_search/uniprot_human_25apr2019.fasta.trypsin.z1-4_nce33.dlib");	
		File smaller=new File("/Users/searleb/Documents/damien/dda_library_search/hela/23aug2017_hela_serum_timecourse_pool_dda_001.dia.dlib");
		//smaller=new File("/Users/searleb/Documents/teaching/encyclopedia/final/quantitative_samples/23aug2017_hela_serum_timecourse_wide_1a.mzML.elib");
		smaller=new File("/Users/searleb/Documents/teaching/encyclopedia/final/example_DIA/23aug2017_hela_serum_timecourse_wide_1a.mzML.elib");

		LibraryFile library=new LibraryFile();
		library.openFile(larger);
		
		LibraryFile file=new LibraryFile();
		file.openFile(smaller);
		ArrayList<LibraryEntry> entries=file.getAllEntries(false, parameters.getAAConstants());
		

		for (LibraryEntry entry : entries) {
			ArrayList<LibraryEntry> candidates=library.getEntries(entry.getPeptideModSeq(), entry.getPrecursorCharge(), false);
			if (candidates.size()>0) {
				AnnotatedLibraryEntry dda=AnnotatedLibraryEntry.getAnnotationsOnly(entry, parameters);
				float correlation=(float)Correlation.getPearsons(candidates.get(0), dda, parameters.getFragmentTolerance());
				
				System.out.println(entry.getPeptideModSeq()+"\t"+entry.getPrecursorCharge()+"\t"+entry.getRetentionTime()/60f+"\t"+candidates.get(0).getRetentionTime()+"\t"+correlation);
				
			}
		}
		
		file.close();
		library.close();
	}
}
