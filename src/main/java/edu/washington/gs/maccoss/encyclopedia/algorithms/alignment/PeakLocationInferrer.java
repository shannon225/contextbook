package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.TreeMap;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.PeptideQuantExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ChromatogramLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakFrequencyCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class PeakLocationInferrer {

	public static PeakLocationInferrerInterface getAlignmentData(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs, ArrayList<PercolatorPeptide> passingPeptides, SearchParameters params) {
		ProgressIndicator subProgress1=new SubProgressIndicator(progress, 0.5f);
		Pair<HashMap<SearchJobData,ArrayList<ChromatogramLibraryEntry>>, HashMap<String,double[]>> pair=getArchetypalPeptides(subProgress1, pecanJobs, passingPeptides, params);
		HashMap<SearchJobData, ArrayList<ChromatogramLibraryEntry>> archetypalPeptides=pair.x;
		HashMap<String, double[]> bestIons=pair.y;

		// get best job
		SearchJobData bestJob=null;
		int max=-1;
		for (Entry<SearchJobData, ArrayList<ChromatogramLibraryEntry>> entry : archetypalPeptides.entrySet()) {
			int length=entry.getValue().size();
			if (length>max) {
				max=length;
				bestJob=entry.getKey();
			}
		}

		Logger.logLine("Seed experiment: "+bestJob.getDiaFile().getName()+" ("+max+" archetypal peptides)");
		Logger.logLine("Seed Percolator file: "+bestJob.getPercolatorFiles().getPeptideOutputFile().getAbsolutePath());
		ArrayList<PercolatorPeptide> alignmentSeed=PercolatorReader.getPassingPeptidesFromTSV(bestJob.getPercolatorFiles().getPeptideOutputFile(), bestJob.getParameters(), false).x;
		TObjectFloatHashMap<String> rtsBySequence=new TObjectFloatHashMap<String>();
		for (PercolatorPeptide peptide : alignmentSeed) {
			rtsBySequence.put(peptide.getPeptideModSeq(), peptide.getRT());
		}
		Logger.logLine("Number of anchors in seed file: "+alignmentSeed.size());

		ProgressIndicator subProgress2=new SubProgressIndicator(progress, 0.5f);

		// construct alignments
		HashMap<SearchJobData, RetentionTimeAlignmentInterface> alignmentMap=new HashMap<SearchJobData, RetentionTimeAlignmentInterface>();
		HashMap<String, Float> alignedRTInMinBySequenceMap=new HashMap<String, Float>();
		// add all bestJob archetypals
		ArrayList<ChromatogramLibraryEntry> archetypals=archetypalPeptides.get(bestJob);
		for (ChromatogramLibraryEntry entry : archetypals) {
			float alignedRT=entry.getRetentionTime()/60f; // no alignment necessary
			alignedRTInMinBySequenceMap.put(entry.getPeptideModSeq(), alignedRT);
		}
		
		int count=0;
		for (SearchJobData job : pecanJobs) {
			if (job!=bestJob) {
				subProgress2.update(job.getDiaFile().getName()+": RT aligning to seed", count/(float) pecanJobs.size());
				count++;

				ArrayList<XYPoint> points=new ArrayList<XYPoint>();
				ArrayList<PercolatorPeptide> peptides=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideOutputFile(), job.getParameters(), false).x;
				for (PercolatorPeptide peptide : peptides) {
					String seq=peptide.getPeptideModSeq();
					if (rtsBySequence.containsKey(seq)) {
						points.add(new XYPoint(rtsBySequence.get(seq)/60f, peptide.getRT()/60f));
					}
				}
				if (points.size()<10) {
					Logger.errorLine("Not enough points ("+points.size()+" out of align:"+peptides.size()+" and best:"+rtsBySequence.size()+") to compute regression between samples, still trying anyways.");
				}
				
				RetentionTimeAlignmentInterface alignment=new RetentionTimeFilter(points, bestJob.getDiaFile().getName(), job.getDiaFile().getName());
				alignmentMap.put(job, alignment);
				if (job instanceof QuantitativeSearchJobData) {
					// try reading encyclopedia data directly from results library
					File resultLibrary=((QuantitativeSearchJobData) job).getResultLibrary();
					alignment.plot(points, Optional.ofNullable(resultLibrary));
				}

				// align local archetypals to the seed
				archetypals=archetypalPeptides.get(job);
				for (ChromatogramLibraryEntry entry : archetypals) {
					float alignedRT=alignment.getXValue(entry.getRetentionTime()/60f);
					alignedRTInMinBySequenceMap.put(entry.getPeptideModSeq(), alignedRT);
				}
			}
		}

		return new SimplePeakLocationInferrer(alignmentMap, new HashMap<>(), alignedRTInMinBySequenceMap, bestIons, params);
	}

	/**
	 * divvy up all the globally passing peptides by the individual search that
	 * best identified them
	 * 
	 * @param progress
	 * @param pecanJobs
	 * @param passingPeptides
	 * @return
	 */
	static Pair<HashMap<SearchJobData, ArrayList<ChromatogramLibraryEntry>>, HashMap<String, double[]>> getArchetypalPeptides(ProgressIndicator progress, ArrayList<SearchJobData> pecanJobs,
			ArrayList<PercolatorPeptide> passingPeptides, SearchParameters params) {
		int numberOfQuantitativePeaks=params.getEffectiveNumberOfQuantitativePeaks();
		MassTolerance fragmentTolerance=params.getFragmentTolerance();
		
		// set up data structures
		HashMap<String, SearchJobData> jobsByFile=new HashMap<String, SearchJobData>();
		HashMap<String, ArrayList<PercolatorPeptide>> peptidesByFile=new HashMap<String, ArrayList<PercolatorPeptide>>();
		for (SearchJobData job : pecanJobs) {
			String name=job.getDiaFile().getName();
			name=name.substring(0, name.lastIndexOf('.'));
			jobsByFile.put(name, job);
			peptidesByFile.put(name, new ArrayList<PercolatorPeptide>());
		}

		// the best individual search is imbedded in the psmID
		for (PercolatorPeptide psm : passingPeptides) {
			String name=psm.getFile();
			name=name.substring(0, name.lastIndexOf('.')); // BRITTLE! Assumes extensions for raw files
			ArrayList<PercolatorPeptide> list=peptidesByFile.get(name);
			if (list==null) {
				Logger.errorLine("Unexpected file ["+name+"] when parsing Percolator result! Ignoring peptide.");
			} else {
				list.add(psm);
			}
		}

		HashMap<String, PeakFrequencyCalculator> ionCounter=new HashMap<String, PeakFrequencyCalculator>();
		
		// extract out chromatogram library entries
		HashMap<SearchJobData, ArrayList<ChromatogramLibraryEntry>> archetypalPeptides=new HashMap<SearchJobData, ArrayList<ChromatogramLibraryEntry>>();
		float increment=1.0f/pecanJobs.size();
		for (Entry<String, ArrayList<PercolatorPeptide>> entry : peptidesByFile.entrySet()) {
			ProgressIndicator subProgress=new SubProgressIndicator(progress, increment);
			
			SearchJobData job=jobsByFile.get(entry.getKey());
			
			ArrayList<PercolatorPeptide> targetPeptides=entry.getValue();

			boolean readFromLibraryResult=false;
			if (job instanceof QuantitativeSearchJobData) {
				// try reading encyclopedia data directly from results library
				File resultLibrary=((QuantitativeSearchJobData) job).getResultLibrary();
				try {
					LibraryInterface results=BlibToLibraryConverter.getFile(resultLibrary);
					ArrayList<LibraryEntry> entries=results.getAllEntries(false, params.getAAConstants());
					ArrayList<ChromatogramLibraryEntry> bestEntries=new ArrayList<ChromatogramLibraryEntry>();

					TreeMap<PeptidePrecursor, LibraryEntry> fastLookupPeptides=new TreeMap<PeptidePrecursor, LibraryEntry>();
					for (LibraryEntry libraryEntry : entries) {
						fastLookupPeptides.put(libraryEntry, libraryEntry);
					}
					
					ArrayList<PercolatorPeptide> missingPeptides=new ArrayList<PercolatorPeptide>();
					for (PercolatorPeptide peptide : targetPeptides) {
						LibraryEntry libEntry=fastLookupPeptides.get(peptide);
						if (!(libEntry instanceof ChromatogramLibraryEntry)) {
							missingPeptides.add(peptide);
							continue;
						}
						// all results files are saved as chromatogram libraries
						ChromatogramLibraryEntry chrom=(ChromatogramLibraryEntry)libEntry;
						String peptideKey=libEntry.getPeptideModSeq()+"+"+libEntry.getPrecursorCharge();
						PeakFrequencyCalculator bestIonsMap=ionCounter.get(peptideKey);
						if (bestIonsMap==null) {
							bestIonsMap=new PeakFrequencyCalculator(fragmentTolerance);
							ionCounter.put(peptideKey, bestIonsMap);
						}
						double[] masses=chrom.getMassArray();
						float[] intensity=chrom.getIntensityArray();
						float[] correlation=chrom.getCorrelationArray();
						for (int i=0; i<correlation.length; i++) {
							if (correlation[i]>=TransitionRefiner.quantitativeCorrelationThreshold) {
								bestIonsMap.increment(masses[i], intensity[i]);
							}
						}
						bestEntries.add(chrom);
					}
					
					if (missingPeptides.size()>0) {
						Logger.logLine("Found "+bestEntries.size()+" archetypal peptides from individual Percolator reports, extracting "+missingPeptides.size()+" additional archetypal peptides from "+job.getDiaFile().getName()+"...");
						ArrayList<ChromatogramLibraryEntry> extracted=extractFromDIA(subProgress, job, missingPeptides, passingPeptides);
						for (ChromatogramLibraryEntry chrom : extracted) {
							String peptideKey=chrom.getPeptideModSeq()+"+"+chrom.getPrecursorCharge();
							PeakFrequencyCalculator bestIonsMap=ionCounter.get(peptideKey);
							if (bestIonsMap==null) {
								bestIonsMap=new PeakFrequencyCalculator(fragmentTolerance);
								ionCounter.put(peptideKey, bestIonsMap);
							}
							double[] masses=chrom.getMassArray();
							float[] intensity=chrom.getIntensityArray();
							float[] correlation=chrom.getCorrelationArray();
							for (int i=0; i<correlation.length; i++) {
								if (correlation[i]>=TransitionRefiner.quantitativeCorrelationThreshold) {
									bestIonsMap.increment(masses[i], intensity[i]);
								}
							}
							bestEntries.add(chrom);
						}
					}
					
					Logger.logLine(resultLibrary.getName()+" produced Parsed:"+targetPeptides.size()+", BEST:"+bestEntries.size());
					
					archetypalPeptides.put(job, bestEntries);
					
					readFromLibraryResult=true;
				} catch (EncyclopediaException e) {
					Logger.errorLine("Parsing error indicates "+job.getPercolatorFiles().getPeptideOutputFile().getName()+" isn't from Encyclopedia:");
					Logger.errorException(e);
					readFromLibraryResult=false;
				} catch (IOException e) {
					throw new EncyclopediaException("Error parsing results library", e);
				} catch (SQLException e) {
					throw new EncyclopediaException("Error parsing results library", e);
				} catch (DataFormatException e) {
					throw new EncyclopediaException("Error parsing results library", e);
				}
			}
			
			// if we can't read data from a library result file (e.g. Pecan), then read directly from the DIA file
			if (!readFromLibraryResult) {
				Logger.logLine("Extracting "+targetPeptides.size()+" Archetypal Peptides from "+job.getDiaFile().getName()+"...");
				subProgress.update(job.getDiaFile().getName()+": Extracting "+targetPeptides.size()+" Archetypal Peptides", 0.00001f);
				ArrayList<ChromatogramLibraryEntry> extracted=extractFromDIA(subProgress, job, targetPeptides, passingPeptides);
				archetypalPeptides.put(job, extracted);
			}
		}
		HashMap<String,double[]> bestIons=new HashMap<String, double[]>();
		for (Entry<String, PeakFrequencyCalculator> entry : ionCounter.entrySet()) {
			String peptideModSeq=entry.getKey();
			double[] ions=entry.getValue().getTopNMasses(numberOfQuantitativePeaks);
			bestIons.put(peptideModSeq, ions);
		}
		return new Pair<HashMap<SearchJobData, ArrayList<ChromatogramLibraryEntry>>, HashMap<String, double[]>>(archetypalPeptides, bestIons);
	}

	/**
	 * Issue 54
	 * This step involves only reading data out of a .dia file (not writing), so we
	 * can operate on the .dia file in place.
	 * 
	 * @param subProgress
	 * @param job
	 * @param targetPeptides
	 * @param passingPeptides
	 * @return
	 */
	private static ArrayList<ChromatogramLibraryEntry> extractFromDIA(ProgressIndicator subProgress, SearchJobData job, ArrayList<PercolatorPeptide> targetPeptides,
			ArrayList<PercolatorPeptide> passingPeptides) {
		StripeFileInterface stripeFile=StripeFileGenerator.getFile(job.getDiaFile(), job.getParameters(), true);

		LibraryInterface library=null;
		if (job instanceof EncyclopediaJobData) {
			library=((EncyclopediaJobData) job).getLibrary();
		}
		ArrayList<IntegratedLibraryEntry> libraryEntries=PeptideQuantExtractor.parseSearchFeatures(subProgress, job, false, passingPeptides, targetPeptides, Optional.ofNullable((PeakLocationInferrerInterface)null), stripeFile, library,
				job.getParameters());
		ArrayList<ChromatogramLibraryEntry> recast=new ArrayList<ChromatogramLibraryEntry>();
		for (IntegratedLibraryEntry e : libraryEntries) {
			recast.add(e);
		}
		stripeFile.close();
		return recast;
	}
}
