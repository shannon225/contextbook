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
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.gui.pecan.PecanJob;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class PecanToBLIB {
	public static void convert(ProgressIndicator progress, ArrayList<PecanJob> pecanJobs, File blibFile) throws IOException, SQLException {
		BlibFile blib=new BlibFile();
		blib.openFile(blibFile);
		blib.dropIndices();
		
		float increment=1.0f/pecanJobs.size();
		for (int i=0; i<pecanJobs.size(); i++) {
			float localComplete=0.0f;
			progress.update("Reading Percolator Results", (i+localComplete)*increment);
			PecanJob job=pecanJobs.get(i);
			File featureFile=job.getFeatureFile();
			File percolatorFile=job.getOutputFile();
			StripeFile stripeFile=new StripeFile();
			stripeFile.openFile(job.getDiaFile());
			PecanScoringFactory taskFactory=job.getTaskFactory();

			ArrayList<ScoredObject<String>> passingPeptides=PercolatorReader.getPassingPeptides(percolatorFile);

			localComplete=0.1f;
			progress.update("Extracting Spectral Data", (i+localComplete)*increment);
			ArrayList<LibraryEntry> libraryEntries=PecanFeatureReader.parsePecanFeatures(featureFile, passingPeptides, stripeFile, taskFactory);

			localComplete=0.5f;
			progress.update("Writing Skyline BLIB", (i+localComplete)*increment);
			// FIXME write BLIB!
		}
	}

}
