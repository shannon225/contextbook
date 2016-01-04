package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanFeatureReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.gui.pecan.PecanJob;
import edu.washington.gs.maccoss.encyclopedia.gui.pecan.PecanPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class PecanToBLIB {
	public static void main(String[] args) {
		File diaFile=new File("/Users/searleb/Documents/school/projects/pecandata/20150708_Ecoli_0911_25x4mzDIA_700_800.dia");
		File fastaFile=new File("/Users/searleb/Documents/school/projects/pecandata/ecoli_20150911_uniprot_sp.fasta");
		File blibFile=new File("/Users/searleb/Documents/school/projects/pecandata/test.blib");
		PecanJob job=PecanPanel.getJob(diaFile, fastaFile, null, SearchParameterParser.parseParameters(SearchParameterParser.getDefaultParameters()));
		ArrayList<PecanJob> jobs=new ArrayList<PecanJob>();
		jobs.add(job);
		convert(new EmptyProgressIndicator(), jobs, blibFile);
	}
	
	public static void convert(ProgressIndicator progress, ArrayList<PecanJob> pecanJobs, File blibFile) {
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
				PecanJob job=pecanJobs.get(i);
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

				localComplete=0.5f;
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
