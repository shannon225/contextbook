package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Date;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchFeatureReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableConcatenator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;

public class SearchToBLIB {
	public static void convert(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, File blibFile, Optional<LibraryFile> libraryFile) {
		ArrayList<File> featureFiles=new ArrayList<File>();
		SearchJobData representativeJob=null;
		for (int i=0; i<pecanJobs.size(); i++) {
			SearchJobData job=pecanJobs.get(i);
			if (!job.hasBeenRun()) {
				continue;
			}
			if (representativeJob==null) {
				representativeJob=job;
			}
			featureFiles.add(job.getFeatureFile());
		}
		
		if (representativeJob==null) return;
		
		File bigFeatureFile=new File(representativeJob.getFeatureFile().getParentFile(), "concatenated_features.txt");
		File bigPercolatorFile=new File(representativeJob.getFeatureFile().getParentFile(), "concatenated_results.txt");

		SearchParameters parameters=representativeJob.getParameters();
		float threshold=parameters.getPercolatorThreshold();
		try {
			ArrayList<ScoredObject<String>> passingPeptides;
			if (bigPercolatorFile.exists()) {
				passingPeptides=PercolatorReader.getPassingPeptidesFromTSV(bigPercolatorFile, threshold);
			} else {
				TableConcatenator.concatenateTables(featureFiles, bigFeatureFile);
				passingPeptides=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorLocation(), bigFeatureFile, bigPercolatorFile, threshold);
			}
			Logger.logLine("Identified "+passingPeptides.size()+" peptides across all files at a "+(threshold*100.0f)+" FDR threshold.");
			convert(progress, pecanJobs, blibFile, libraryFile, Optional.of(passingPeptides));
			progress.update(passingPeptides.size()+" peptides identified at "+(threshold*100.0f)+"% FDR", 1.0f);
		} catch (IOException ioe) {
			Logger.errorLine("Error creating concatenated feature file");
			Logger.errorException(ioe);
		} catch (InterruptedException ie) {
			Logger.errorLine("Error creating concatenated feature file");
			Logger.errorException(ie);
		}
	}
	
	static void convert(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, File blibFile, Optional<LibraryFile> libraryFile, Optional<ArrayList<ScoredObject<String>>> passingPeptides) {
		try {
			BlibFile blib=new BlibFile();
			blib.openFile();
			blib.setUserFile(blibFile);
			blib.dropIndices();
			int[] counterTotals=new int[] {0,0,0};

			float increment=1.0f/pecanJobs.size();
			for (int i=0; i<pecanJobs.size(); i++) {
				SearchJobData job=pecanJobs.get(i);
				if (!job.hasBeenRun()) {
					continue;
				}
				ProgressIndicator subProgress=new SubProgressIndicator(progress, increment);
				
				ArrayList<ScoredObject<String>> globalPassingPeptides;
				ArrayList<ScoredObject<String>>	localPassingPeptides=PercolatorReader.getPassingPeptidesFromTSV(job.getOutputFile(), pecanJobs.get(i).getParameters().getPercolatorThreshold());
				if (passingPeptides.isPresent()) {
					globalPassingPeptides=passingPeptides.get();
				} else {
					globalPassingPeptides=localPassingPeptides;
				}
				
				counterTotals=convertFile(subProgress, job, globalPassingPeptides, localPassingPeptides, counterTotals, blib, libraryFile);
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

	static int[] convertFile(ProgressIndicator subProgress, SearchJobData job, ArrayList<ScoredObject<String>> globalPassingPeptides, ArrayList<ScoredObject<String>> localPassingPeptides, int[] counterTotals, BlibFile blib, Optional<LibraryFile> libraryFile) throws IOException, SQLException {
		File diaFile=job.getDiaFile();
		Logger.logLine("Reading Percolator Results from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Reading Percolator Results", 0.0f);

		File featureFile=job.getFeatureFile();

		StripeFileInterface stripeFile=MzmlToDIAConverter.getFile(diaFile, job.getParameters());
		Logger.logLine("Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides", 0.1f);

		ArrayList<IntegratedLibraryEntry> libraryEntries=SearchFeatureReader.parseSearchFeatures(featureFile, globalPassingPeptides, localPassingPeptides, stripeFile, libraryFile, job.getParameters());

		float totalTIC=0.0f;
		for (IntegratedLibraryEntry entry : libraryEntries) {
			totalTIC+=entry.getTIC();
		}
		float normalizer=totalTIC/1e12f;
		
		File integrationFile=new File(diaFile.getAbsolutePath()+".integration.txt");

		PrintWriter writer=new PrintWriter(integrationFile, "UTF-8");
		writer.println("File\tPeptideModSeq\tPrecursorCharge\tFragmentIons\tRTStart\tRTCenter\tRTStop\tTIC\tNormTIC");
		for (IntegratedLibraryEntry entry : libraryEntries) {
			writer.println(diaFile.getName()+"\t"+entry.getPeptideModSeq()+"\t"+entry.getPrecursorCharge()+"\t"+entry.getIonCount()+"\t"+entry.getRtRange().getStart()+"\t"+entry.getRetentionTime()+"\t"+entry.getRtRange().getStop()+"\t"+entry.getTIC()+"\t"+(entry.getTIC()/normalizer));
		}
		writer.flush();
		writer.close();
		
		ArrayList<LibraryEntry> recasted=new ArrayList<LibraryEntry>();
		for (IntegratedLibraryEntry entry : libraryEntries) {
			recasted.add(entry);
		}

		Logger.logLine("Writing Skyline BLIB from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Writing Skyline BLIB", 0.9f);

		counterTotals=blib.addLibrary(job, recasted, counterTotals[0], counterTotals[1], counterTotals[2]);
		subProgress.update(diaFile.getName()+": Finished writing to Skyline BLIB at"+new Date().toString(), 1.0f);
		return counterTotals;
	}
}
