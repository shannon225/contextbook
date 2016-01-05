package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;

import com.google.common.base.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanFeatureReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class PecanToBLIB {
	public static void main(String[] args) {
		SearchParameters parameters=SearchParameterParser.parseParameters(SearchParameterParser.getDefaultParameters());
		File fastaFile=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/ecoli-190209-contam_correctNL.fasta");
		File blibFile=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/test.blib");
		ArrayList<FastaEntry> targets=null;
		
		File diaFile1=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/20150708_Ecoli_0911_25x4mzDIA_500_600.dia");
		File outputFile1=new File(diaFile1.getAbsolutePath()+".pecan.txt");
		File featureFile1=new File(outputFile1.getAbsolutePath()+".features.txt");
		PecanScoringFactory factory1=new PecanOneScoringFactory(parameters, featureFile1);
		PecanJobData job1=new PecanJobData(Optional.fromNullable(targets), diaFile1, fastaFile,featureFile1, outputFile1, factory1);

		File diaFile2=new File("/Users/searleb/Documents/projects/pecan/ecoli_dataset/20150708_Ecoli_0911_25x4mzDIA_600_700.dia");
		File outputFile2=new File(diaFile2.getAbsolutePath()+".pecan.txt");
		File featureFile2=new File(outputFile2.getAbsolutePath()+".features.txt");
		PecanScoringFactory factory2=new PecanOneScoringFactory(parameters, featureFile2);
		PecanJobData job2=new PecanJobData(Optional.fromNullable(targets), diaFile2, fastaFile,featureFile2, outputFile2, factory2);
		
		ArrayList<PecanJobData> jobs=new ArrayList<PecanJobData>();
		jobs.add(job1);
		jobs.add(job2);
		
		convert(new EmptyProgressIndicator(), jobs, blibFile);
	}
	
	public static void convert(ProgressIndicator progress, ArrayList<PecanJobData> pecanJobs, File blibFile) {
		try {
			BlibFile blib=new BlibFile();
			blib.openFile();
			blib.setUserFile(blibFile);
			blib.dropIndices();
			int idCounter=0;
			int jobCounter=0;
			int modCounter=0;

			float increment=1.0f/pecanJobs.size();
			for (int i=0; i<pecanJobs.size(); i++) {
				float localComplete=0.0f;
				PecanJobData job=pecanJobs.get(i);
				progress.update(job.getDiaFile().getName()+": Reading Percolator Results", (i+localComplete)*increment);

				File featureFile=job.getFeatureFile();
				File percolatorFile=job.getOutputFile();
				StripeFile stripeFile=new StripeFile();
				stripeFile.openFile(job.getDiaFile());
				PecanScoringFactory taskFactory=job.getTaskFactory();

				ArrayList<ScoredObject<String>> passingPeptides=PercolatorReader.getPassingPeptides(percolatorFile, taskFactory.getParameters().getPercolatorThreshold());

				localComplete=0.1f;
				progress.update(job.getDiaFile().getName()+": Extracting Spectral Data for "+passingPeptides.size()+" Peptides", (i+localComplete)*increment);
				ArrayList<LibraryEntry> libraryEntries=PecanFeatureReader.parsePecanFeatures(featureFile, passingPeptides, stripeFile, taskFactory);

				localComplete=0.9f;
				progress.update(job.getDiaFile().getName()+": Writing Skyline BLIB", (i+localComplete)*increment);

				int[] counters=blib.addLibrary(job, libraryEntries, idCounter, jobCounter, modCounter);
				idCounter=counters[0];
				jobCounter=counters[1];
				modCounter=counters[2];
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

}
