package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Date;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideQuantExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableConcatenator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;

public class SearchToBLIB {
	public static void convert(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, File libFile, boolean writeBlib, Optional<LibraryInterface> libraryFile) {
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
		float threshold=parameters.getEffectivePercolatorThreshold();
		try {
			ArrayList<ScoredObject<String>> passingPeptides;
			if (featureFiles.size()==1) {
				// if there's only one file then don't need to re-run percolator
				passingPeptides=PercolatorReader.getPassingPeptidesFromTSV(representativeJob.getOutputFile(), threshold);
			} else {
				TableConcatenator.concatenateTables(featureFiles, bigFeatureFile);
				passingPeptides=PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorLocation(), bigFeatureFile, bigPercolatorFile, threshold);
			}
			Logger.logLine("Identified "+passingPeptides.size()+" peptides across all files at a "+(threshold*100.0f)+" FDR threshold.");
			if (writeBlib) {
				convertBlib(progress, pecanJobs, libFile, libraryFile, Optional.of(passingPeptides));
			} else {
				convertElib(progress, pecanJobs, libFile, libraryFile, Optional.of(passingPeptides));
			}
			progress.update(passingPeptides.size()+" peptides identified at "+(threshold*100.0f)+"% FDR", 1.0f);
		} catch (IOException ioe) {
			Logger.errorLine("Error creating concatenated feature file");
			Logger.errorException(ioe);
		} catch (InterruptedException ie) {
			Logger.errorLine("Error creating concatenated feature file");
			Logger.errorException(ie);
		}
	}
	
	static void convertBlib(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, File blibFile, Optional<LibraryInterface> libraryFile, Optional<ArrayList<ScoredObject<String>>> passingPeptides) {
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
				ArrayList<ScoredObject<String>>	localPassingPeptides=PercolatorReader.getPassingPeptidesFromTSV(job.getOutputFile(), pecanJobs.get(i).getParameters().getEffectivePercolatorThreshold());
				if (passingPeptides.isPresent()) {
					globalPassingPeptides=passingPeptides.get();
				} else {
					globalPassingPeptides=localPassingPeptides;
				}
				
				counterTotals=convertFileBlib(subProgress, job, globalPassingPeptides, localPassingPeptides, counterTotals, blib, libraryFile);
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

	/**
	 * trims to quantifiable peptides! for loading into skyline!
	 * @param subProgress
	 * @param job
	 * @param globalPassingPeptides
	 * @param localPassingPeptides
	 * @param counterTotals
	 * @param blib
	 * @param libraryFile
	 * @return
	 * @throws IOException
	 * @throws SQLException
	 */
	static int[] convertFileBlib(ProgressIndicator subProgress, SearchJobData job, ArrayList<ScoredObject<String>> globalPassingPeptides, ArrayList<ScoredObject<String>> localPassingPeptides, int[] counterTotals, BlibFile blib, Optional<LibraryInterface> libraryFile) throws IOException, SQLException {
		File diaFile=job.getDiaFile();
		Logger.logLine("Reading Percolator Results from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Reading Percolator Results", 0.0f);

		File featureFile=job.getFeatureFile();

		StripeFileInterface stripeFile=MzmlToDIAConverter.getFile(diaFile, job.getParameters());
		Logger.logLine("Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides", 0.00001f);

		//ArrayList<IntegratedLibraryEntry> libraryEntries=SearchFeatureReader.parseSearchFeatures(featureFile, globalPassingPeptides, localPassingPeptides, stripeFile, Optional.ofNullable((LibraryFile)null), job.getParameters());
		ArrayList<IntegratedLibraryEntry> libraryEntries=PeptideQuantExtractor.parseSearchFeatures(subProgress, featureFile, true, globalPassingPeptides, localPassingPeptides, stripeFile, libraryFile, job.getParameters());
		stripeFile.close();
		
		File integrationFile=new File(diaFile.getAbsolutePath()+".integration.txt");

		PrintWriter writer=new PrintWriter(integrationFile, "UTF-8");
		writer.println("File Name\tPeptide Modified Sequence\tMin Start Time\tMax End Time\tPrecursor Charge\tPrecursorIsDecoy\tIon Count\tRetention Time Center\tTIC");
		
		for (IntegratedLibraryEntry entry : libraryEntries) {
			writer.println(diaFile.getName()+"\t"+entry.getPeptideModSeq()+"\t"+entry.getRtRange().getStart()/60f+"\t"+entry.getRtRange().getStop()/60f+"\t"+entry.getPrecursorCharge()+"\tFALSE\t"+entry.getIonCount()+"\t"+entry.getRetentionTime()/60f+"\t"+entry.getTIC());
		}
		writer.flush();
		writer.close();
		
		ArrayList<LibraryEntry> recasted=new ArrayList<LibraryEntry>();
		for (IntegratedLibraryEntry entry : libraryEntries) {
			recasted.add(entry);
		}

		Logger.logLine("Writing Skyline BLIB from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Writing Skyline BLIB", 0.99999f);

		counterTotals=blib.addLibrary(job, recasted, counterTotals[0], counterTotals[1], counterTotals[2]);
		subProgress.update(diaFile.getName()+": Finished writing to Skyline BLIB at"+new Date().toString(), 1.0f);
		return counterTotals;
	}
	
	static void convertElib(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, File elibFile, Optional<LibraryInterface> libraryFile, Optional<ArrayList<ScoredObject<String>>> passingPeptides) {
		try {
			LibraryFile elib=new LibraryFile();
			elib.openFile();
			elib.dropIndices();

			float increment=1.0f/pecanJobs.size();
			for (int i=0; i<pecanJobs.size(); i++) {
				SearchJobData job=pecanJobs.get(i);
				if (!job.hasBeenRun()) {
					continue;
				}
				ProgressIndicator subProgress=new SubProgressIndicator(progress, increment);
				
				ArrayList<ScoredObject<String>> globalPassingPeptides;
				ArrayList<ScoredObject<String>>	localPassingPeptides=PercolatorReader.getPassingPeptidesFromTSV(job.getOutputFile(), pecanJobs.get(i).getParameters().getEffectivePercolatorThreshold());
				if (passingPeptides.isPresent()) {
					globalPassingPeptides=passingPeptides.get();
				} else {
					globalPassingPeptides=localPassingPeptides;
				}
				
				convertFileElib(subProgress, job, globalPassingPeptides, localPassingPeptides, elib, libraryFile);
			}

			elib.createIndices();
			elib.saveAsFile(elibFile);
			
		} catch (IOException ioe) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(ioe);
		} catch (SQLException sqle) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(sqle);
		}
	}

	/**
	 * Does not limit to quantifiable! Reports all potential peaks!
	 * @param subProgress
	 * @param job
	 * @param globalPassingPeptides
	 * @param localPassingPeptides
	 * @param elib
	 * @param libraryFile
	 * @throws IOException
	 * @throws SQLException
	 */
	static void convertFileElib(ProgressIndicator subProgress, SearchJobData job, ArrayList<ScoredObject<String>> globalPassingPeptides, ArrayList<ScoredObject<String>> localPassingPeptides, LibraryFile elib, Optional<LibraryInterface> libraryFile) throws IOException, SQLException {
		File diaFile=job.getDiaFile();
		Logger.logLine("Reading Percolator Results from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Reading Percolator Results", 0.0f);

		File featureFile=job.getFeatureFile();

		StripeFileInterface stripeFile=MzmlToDIAConverter.getFile(diaFile, job.getParameters());
		Logger.logLine("Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides", 0.00001f);

		ArrayList<IntegratedLibraryEntry> libraryEntries=PeptideQuantExtractor.parseSearchFeatures(subProgress, featureFile, false, globalPassingPeptides, localPassingPeptides, stripeFile, libraryFile, job.getParameters());
		stripeFile.close();
		
		Logger.logLine("Writing Encyclopedia ELIB from "+diaFile.getName()+"...");
		subProgress.update(diaFile.getName()+": Writing Encyclopedia ELIB", 0.99999f);

		elib.addIntegratedEntries(libraryEntries, TransitionRefiner.identificationCorrelationThreshold);
		subProgress.update(diaFile.getName()+": Finished writing to Encyclopedia ELIB at"+new Date().toString(), 1.0f);
	}
}
