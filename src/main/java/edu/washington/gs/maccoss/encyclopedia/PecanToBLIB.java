package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Date;

import com.google.common.base.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanFeatureReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableConcatenator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;

public class PecanToBLIB {
	public static void main(String[] args) {
		PecanSearchParameters parameters=PecanParameterParser.parseParameters(PecanParameterParser.getDefaultParameters());
		File fastaFile=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/ecoli-190209-contam_correctNL.fasta");
		File blibFile=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/ecoli.blib");
		ArrayList<FastaEntry> targets=null;
		
		PecanJobData job1=getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/pecan/ecoli_dataset/20150708_Ecoli_0911_25x4mzDIA_500_600.dia");
		PecanJobData job2=getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/pecan/ecoli_dataset/20150708_Ecoli_0911_25x4mzDIA_600_700.dia");
		PecanJobData job3=getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/pecan/ecoli_dataset/20150708_Ecoli_0911_25x4mzDIA_700_800.dia");
		PecanJobData job4=getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/pecan/ecoli_dataset/20150708_Ecoli_0911_25x4mzDIA_800_900.dia");

		ArrayList<PecanJobData> jobs=new ArrayList<PecanJobData>();
		jobs.add(job1);
		jobs.add(job2);
		jobs.add(job3);
		jobs.add(job4);
		
		convert(new EmptyProgressIndicator(), jobs, blibFile);
	}

	static PecanJobData getData(PecanSearchParameters parameters, File fastaFile, ArrayList<FastaEntry> targets, String dia) {
		File diaFile1=new File(dia);
		File outputFile1=new File(diaFile1.getAbsolutePath()+".pecan.txt");
		File featureFile1=new File(outputFile1.getAbsolutePath()+".features.txt");
		PecanScoringFactory factory1=new PecanOneScoringFactory(parameters, featureFile1);
		PecanJobData job1=new PecanJobData(Optional.fromNullable(targets), diaFile1, fastaFile,featureFile1, outputFile1, factory1);
		return job1;
	}

	public static void convert(ProgressIndicator progress, ArrayList<PecanJobData> pecanJobs, File blibFile) {
		ArrayList<File> featureFiles=new ArrayList<File>();
		PecanJobData representativeJob=null;
		for (int i=0; i<pecanJobs.size(); i++) {
			PecanJobData job=pecanJobs.get(i);
			if (!job.hasBeenRun()) {
				continue;
			}
			if (representativeJob==null) {
				representativeJob=job;
			}
			featureFiles.add(job.getFeatureFile());
		}
		
		if (representativeJob==null) return;
		
		File bigFeatureFile=new File(representativeJob.getFeatureFile().getParentFile(), "concatenated_pecan_features.txt");
		File bigPercolatorFile=new File(representativeJob.getFeatureFile().getParentFile(), "concatenated_pecan_results.txt");
		
		try {
			TableConcatenator.concatenateTables(featureFiles, bigFeatureFile);
			SearchParameters parameters=representativeJob.getTaskFactory().getParameters();
			float threshold=parameters.getPercolatorThreshold();
			ArrayList<ScoredObject<String>> passingPeptides=PercolatorExecutor.executePercolator(parameters.getPercolatorLocation(), bigFeatureFile, bigPercolatorFile, threshold);
			Logger.logLine("Identified "+passingPeptides.size()+" peptides across all files at a "+(threshold*100.0f)+" FDR threshold.");
			convert(progress, pecanJobs, blibFile, Optional.of(passingPeptides));
			progress.update(passingPeptides.size()+" peptides identified at "+(threshold*100.0f)+"% FDR", 1.0f);
		} catch (IOException ioe) {
			Logger.errorLine("Error creating concatenated feature file");
			Logger.errorException(ioe);
		} catch (InterruptedException ie) {
			Logger.errorLine("Error creating concatenated feature file");
			Logger.errorException(ie);
		}
	}
	
	static void convert(ProgressIndicator progress, ArrayList<PecanJobData> pecanJobs, File blibFile, Optional<ArrayList<ScoredObject<String>>> passingPeptides) {
		try {
			BlibFile blib=new BlibFile();
			blib.openFile();
			blib.setUserFile(blibFile);
			blib.dropIndices();
			int[] counterTotals=new int[] {0,0,0};

			float increment=1.0f/pecanJobs.size();
			for (int i=0; i<pecanJobs.size(); i++) {
				PecanJobData job=pecanJobs.get(i);
				if (!job.hasBeenRun()) {
					continue;
				}
				ProgressIndicator subProgress=new SubProgressIndicator(progress, increment);
				
				ArrayList<ScoredObject<String>> localPassingPeptides;
				if (passingPeptides.isPresent()) {
					localPassingPeptides=passingPeptides.get();
				} else {
					localPassingPeptides=PercolatorReader.getPassingPeptides(job.getOutputFile(), pecanJobs.get(i).getTaskFactory().getParameters().getPercolatorThreshold());
				}
				
				counterTotals=convertFile(subProgress, job, localPassingPeptides, counterTotals, blib);
			}

			blib.createIndices();
			blib.saveFile();
		} catch (IOException ioe) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(ioe);
		} catch (SQLException sqle) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(sqle);
		}
	}

	static int[] convertFile(ProgressIndicator subProgress, PecanJobData job, ArrayList<ScoredObject<String>> passingPeptides, int[] counterTotals, BlibFile blib) throws IOException, SQLException {
		File diaFile=job.getDiaFile();
		Logger.logLine("Reading Percolator Results from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Reading Percolator Results", 0.0f);

		File featureFile=job.getFeatureFile();

		PecanScoringFactory taskFactory=job.getTaskFactory();
		StripeFileInterface stripeFile=MzmlToDIAConverter.getFile(diaFile, taskFactory.getParameters());
		Logger.logLine("Extracting Spectral Data for "+passingPeptides.size()+" Peptides from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Extracting Spectral Data for "+passingPeptides.size()+" Peptides", 0.1f);

		ArrayList<LibraryEntry> libraryEntries=PecanFeatureReader.parsePecanFeatures(featureFile, passingPeptides, stripeFile, taskFactory);

		Logger.logLine("Writing Skyline BLIB from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Writing Skyline BLIB", 0.9f);

		counterTotals=blib.addLibrary(job, libraryEntries, counterTotals[0], counterTotals[1], counterTotals[2]);
		subProgress.update(diaFile.getName()+": Finished writing to Skyline BLIB at"+new Date().toString(), 1.0f);
		return counterTotals;
	}
}
