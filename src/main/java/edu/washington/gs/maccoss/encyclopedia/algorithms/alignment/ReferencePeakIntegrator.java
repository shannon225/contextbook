package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ParsimonyProteinGrouper;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface.AlignmentDataPoint;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.PeptideQuantExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.PeptideQuantExtractorTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.RelativePeakIntensityMatrix;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ChromatogramLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DDASearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.LimitedQueue;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class ReferencePeakIntegrator {
	private final SearchParameters params;
	private final HashMap<SearchJobData, RetentionTimeAlignmentInterface> alignmentsBySample;
	private final HashMap<SearchJobData, List<AlignmentDataPoint>> pointsBySample;
	private final HashMap<String, LibraryEntry> peaksByPeptide;
	private final HashMap<String, PSMData> psmDataByPeptide;
	
	public ReferencePeakIntegrator(HashMap<SearchJobData, RetentionTimeAlignmentInterface> alignmentsBySample, HashMap<SearchJobData, List<AlignmentDataPoint>> pointsBySample,
			HashMap<String, LibraryEntry> peaksByPeptide, HashMap<String, PSMData> psmDataByPeptide, SearchParameters params) {
		this.alignmentsBySample = alignmentsBySample;
		this.pointsBySample = pointsBySample;
		this.peaksByPeptide = peaksByPeptide;
		this.psmDataByPeptide=psmDataByPeptide;
		this.params=params;
	}
	
	public static LibraryFile integrateAllPeptides(File fileToWrite, LibraryInterface reference, List<? extends SearchJobData> jobs, ArrayList<PercolatorPeptide> passingPeptides, SearchParameters params, ProgressIndicator progress) {
		ProgressIndicator subProgress=new SubProgressIndicator(progress, 0.1f);
		ReferencePeakIntegrator integrator;

		try {
			integrator=buildIntegrator(reference, jobs, passingPeptides, params, subProgress);
		} catch (IOException ioe) {
			throw new EncyclopediaException("error building integrator", ioe);
		} catch (SQLException sqle) {
			throw new EncyclopediaException("error building integrator", sqle);
		} catch (DataFormatException dfe) {
			throw new EncyclopediaException("error building integrator", dfe);
		}
		
		HashMap<String, RelativePeakIntensityMatrix> matricies=new HashMap<>();
		

		try {
			LibraryFile elib=new LibraryFile();
			elib.openFile();
			elib.dropIndices();
	
			ArrayList<String> keys=new ArrayList<>();
			for (SearchJobData job : jobs) {
				String sampleName=job.getDiaFileReader().getOriginalFileName();
				sampleName=sampleName.substring(0, sampleName.lastIndexOf('.')); // file names are lose extensions
				keys.add(sampleName);
				try {
					ArrayList<IntegratedLibraryEntry> results=integrator.integratePeptides(job);
					elib.addIntegratedEntries(!(job instanceof DDASearchJobData), results, Optional.empty(), Optional.empty(), params);
					for (IntegratedLibraryEntry entry : results) {
						RelativePeakIntensityMatrix matrix=matricies.get(entry.getPeptideModSeq());
						if (matrix==null) {
							LibraryEntry target=integrator.peaksByPeptide.get(entry.getPeptideModSeq());
							matrix=new RelativePeakIntensityMatrix(params.getFragmentTolerance(), Optional.of(target));
							matricies.put(entry.getPeptideModSeq(), matrix);
						}
						matrix.addPeak(sampleName, entry);
					}
					
					progress.update("Finished processing "+job.getDiaFileReader().getOriginalFileName(), 0.9f/jobs.size());
					Logger.logLine("Finished processing "+job.getDiaFileReader().getOriginalFileName()+", found "+results.size()+" peptides.");
					
				} catch (InterruptedException ie) {
					throw new EncyclopediaException("error using integrator on ["+job.getDiaFileReader().getOriginalFileName()+"]", ie);
				} catch (IOException ioe) {
					throw new EncyclopediaException("error using integrator on ["+job.getDiaFileReader().getOriginalFileName()+"]", ioe);
				} catch (SQLException sqle) {
					throw new EncyclopediaException("error using integrator on ["+job.getDiaFileReader().getOriginalFileName()+"]", sqle);
				} catch (DataFormatException dfe) {
					throw new EncyclopediaException("error using integrator on ["+job.getDiaFileReader().getOriginalFileName()+"]", dfe);
				}
			}
			
			elib.addProteinsFromPercolator(passingPeptides);
			elib.addMetadata(params.toParameterMap());
			elib.setSources(jobs);

			elib.createIndices();
			elib.saveAsFile(fileToWrite);
		
			String saveFilePrefix=fileToWrite.getAbsolutePath();
			File f=new File(saveFilePrefix+".pepquant.txt");
			PrintWriter writer=new PrintWriter(f, "UTF-8");
			
			Collections.sort(keys);
			StringBuilder sb=new StringBuilder("peptide\tnumPeaks");
			for (String key : keys) {
				sb.append("\t");
				sb.append(key);
			}
			writer.println(sb);
			
			ArrayList<String> sortedPeptides=new ArrayList<>(matricies.keySet());
			Collections.sort(sortedPeptides);
			for (String peptideModSeq : sortedPeptides) {
				RelativePeakIntensityMatrix matrix=matricies.get(peptideModSeq);
				double[] nBestPeaks=matrix.pickNBestPeaks(params.getNumberOfQuantitativePeaks(), params.getMinNumOfQuantitativePeaks());
	
				if (nBestPeaks.length>=params.getMinNumOfQuantitativePeaks()&&nBestPeaks.length>0) {
					float[] intensities=matrix.integratePeptides(nBestPeaks, keys);
	
					sb.setLength(0);
					sb.append(peptideModSeq);
					sb.append("\t");
					sb.append(nBestPeaks.length);
					
					for (int i = 0; i < intensities.length; i++) {
						sb.append("\t");
						sb.append(intensities[i]);
					}
					writer.println(sb);
				}
			}
			writer.close();
			return elib;
		} catch (IOException ioe) {
			throw new EncyclopediaException("error writing integrator data!", ioe);
		} catch (SQLException sqle) {
			throw new EncyclopediaException("error writing integrator data!", sqle);
		}
	}
	
	public ArrayList<IntegratedLibraryEntry> integratePeptides(SearchJobData job) throws IOException, SQLException, DataFormatException, InterruptedException {
		ConcurrentLinkedQueue<IntegratedLibraryEntry> savedEntries=new ConcurrentLinkedQueue<IntegratedLibraryEntry>();
		RetentionTimeAlignmentInterface alignment=alignmentsBySample.get(job);
		
		StripeFileInterface stripeFile = job.getDiaFileReader();
		String filename=stripeFile.getOriginalFileName();
		

		ArrayList<Range> ranges=new ArrayList<Range>();
		for (Range range : stripeFile.getRanges().keySet()) {
			ranges.add(range);
		}
		Collections.sort(ranges);
		
		ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("REFERENCE_PEAK_INTEGRATOR-%d").setDaemon(true).build();
		LimitedQueue<Runnable> workQueue=new LimitedQueue<Runnable>(10000);
		ExecutorService executor=new ThreadPoolExecutor(params.getNumberOfThreadsUsed(), params.getNumberOfThreadsUsed(), Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

		for (Range range : ranges) {
			boolean used=false;
			float minRetentionTime=Float.MAX_VALUE;
			float maxRetentionTime=-Float.MAX_VALUE;
			Collection<PSMData> targets = psmDataByPeptide.values();
			for (PSMData psm : targets) {
				if (range.contains((float)psm.getPrecursorMZ())) {
					minRetentionTime=Math.min(minRetentionTime, psm.getRetentionTime()-psm.getDuration());
					maxRetentionTime=Math.max(maxRetentionTime, psm.getRetentionTime()+psm.getDuration());
					used=true;
				}
			}
			if (!used) continue;
			String baseMessage="Extracting "+range+" m/z ("+Math.max(0.0f, minRetentionTime/60f)+" to "+Math.max(0.0f, maxRetentionTime/60f)+" min)";
			Logger.logLine("Quant "+baseMessage);

			ArrayList<FragmentScan> stripes=stripeFile.getStripes(range.getMiddle(), minRetentionTime, maxRetentionTime, false);
			Collections.sort(stripes);

			// prepare executor for background
			for (PSMData psm : targets) {
				if (range.contains((float)psm.getPrecursorMZ())) {
					// convert RTs back to this sample
					PSMData rtCorrected=psm.updateRetentionTime(alignment.getYValue(psm.getRetentionTime()));
					executor.submit(new PeptideQuantExtractorTask(filename, rtCorrected, Optional.empty(), Optional.empty(), stripes, params, savedEntries, false));
				}
			}
		}

		executor.shutdown();
		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

		ArrayList<IntegratedLibraryEntry> entryList=new ArrayList<IntegratedLibraryEntry>();
		for (IntegratedLibraryEntry entry : savedEntries) {
			entryList.add(entry);
		}

		return entryList;
	}

	public static ReferencePeakIntegrator buildIntegrator(LibraryInterface reference, List<? extends SearchJobData> jobs, ArrayList<PercolatorPeptide> passingPeptides, SearchParameters params, ProgressIndicator progress) throws IOException, SQLException, DataFormatException {
		HashMap<SearchJobData, RetentionTimeAlignmentInterface> alignmentsBySample=new HashMap<>();
		HashMap<SearchJobData, List<AlignmentDataPoint>> pointsBySample=new HashMap<>();
		HashMap<String, LibraryEntry> targetByPeptide=new HashMap<>();
		HashMap<String, PSMData> psmDataByPeptide=new HashMap<>();

		HashMap<String, SearchJobData> jobsByFile=new HashMap<>();
		HashMap<SearchJobData, HashMap<String, PercolatorPeptide>> peptidesByFile=new HashMap<>();
		for (SearchJobData job : jobs) {
			String name=job.getDiaFileReader().getOriginalFileName();
			name=name.substring(0, name.lastIndexOf('.')); // file names are lose extensions
			jobsByFile.put(name, job);
			peptidesByFile.put(job, new HashMap<String, PercolatorPeptide>());
		}
		
		// the best individual search is imbedded in the psmID
		for (PercolatorPeptide psm : passingPeptides) {
			String name=psm.getFile();
			name=name.substring(0, name.lastIndexOf('.')); // file names are lose extensions
			SearchJobData key = jobsByFile.get(name);
			if (key==null) {
				Logger.errorLine("Unexpected file ["+name+"] when parsing Percolator result! Ignoring peptide.");
			} else {
				HashMap<String, PercolatorPeptide> list=peptidesByFile.get(key);
				list.put(psm.getPeptideModSeq(), psm);
			}
		}
		
		TObjectFloatHashMap<String> referenceRTInSecs=getRTs(reference, params);

		for (SearchJobData job : jobs) {
			progress.update("Aligning "+job.getDiaFileReader().getOriginalFileName()+" to reference...");
			
			if (job instanceof QuantitativeSearchJobData) {
				HashMap<String, PercolatorPeptide> thisSamplesTargets=peptidesByFile.get(job);
				File resultLibrary=((QuantitativeSearchJobData) job).getResultLibrary();
				LibraryInterface results=BlibToLibraryConverter.getFile(resultLibrary);
				
				Pair<RetentionTimeAlignmentInterface, List<AlignmentDataPoint>> pair=getAlignment(results, params, referenceRTInSecs, resultLibrary);
				alignmentsBySample.put(job, pair.x);
				pointsBySample.put(job, pair.y);

				// pull targets from the elib if possible (RTs are converted to reference)
				int count=0;
				for (LibraryEntry entry : results.getAllEntries(false, params.getAAConstants())) {
					PercolatorPeptide target=thisSamplesTargets.get(entry.getPeptideModSeq());
					if (target!=null) {
						LibraryEntry rtCorrectedEntry=entry.updateRetentionTime(pair.x.getXValue(entry.getRetentionTime()));
						targetByPeptide.put(rtCorrectedEntry.getPeptideModSeq(), rtCorrectedEntry);
						
						thisSamplesTargets.remove(entry.getPeptideModSeq());
						count++;
					}
				}	
				Logger.logLine("Targets in "+resultLibrary.getName()+": found "+count+" in individual search, still need to look up "+thisSamplesTargets.size());
				
				// otherwise we need to look them up in the file (RTs are converted to reference)
				if (thisSamplesTargets.size()>0) {
					ArrayList<PercolatorPeptide> localPassingPeptides = new ArrayList<>(thisSamplesTargets.values());
					LibraryInterface library=null;
					if (job instanceof EncyclopediaJobData) {
						library=((EncyclopediaJobData)job).getLibrary();
					}
					StripeFileInterface stripeFile = job.getDiaFileReader();
					ArrayList<IntegratedLibraryEntry> libraryEntries=PeptideQuantExtractor.parseSearchFeatures(new EmptyProgressIndicator(false), job, false, localPassingPeptides, localPassingPeptides, Optional.empty(), stripeFile, library, job.getParameters());
					for (IntegratedLibraryEntry entry : libraryEntries) {
						LibraryEntry rtCorrectedEntry=entry.updateRetentionTime(pair.x.getXValue(entry.getRetentionTime()));
						targetByPeptide.put(rtCorrectedEntry.getPeptideModSeq(), rtCorrectedEntry);
						thisSamplesTargets.remove(entry.getPeptideModSeq());
					}
				}
				if (thisSamplesTargets.size()>0) {
					Logger.errorLine("Targets in "+resultLibrary.getName()+" remaining: "+thisSamplesTargets.size());
				}
				
			}
		}
		
		for (LibraryEntry entry : targetByPeptide.values()) {
			float duration;
			if (entry instanceof ChromatogramLibraryEntry) {
				duration=((ChromatogramLibraryEntry) entry).getDurationInSec();
			} else {
				duration=params.getExpectedPeakWidth();
			}
			PSMData psmdata=new PSMData(entry.getAccessions(), entry.getSpectrumIndex(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), entry.getPeptideModSeq(), entry.getRetentionTime(), entry.getScore(), 1.0f-entry.getScore(), duration, false, params.getAAConstants());
			psmDataByPeptide.put(entry.getPeptideModSeq(), psmdata);
		}
		
		Logger.logLine("Found "+targetByPeptide.size()+" total reference peptides.");
		
		return new ReferencePeakIntegrator(alignmentsBySample, pointsBySample, targetByPeptide, psmDataByPeptide, params);
	}

	private static Pair<RetentionTimeAlignmentInterface, List<AlignmentDataPoint>> getAlignment(LibraryInterface library, SearchParameters params, TObjectFloatHashMap<String> referenceRTInSecs, File saveFileSeed) throws IOException, SQLException, DataFormatException {
		TObjectFloatHashMap<String> wideRTInSecs = getRTs(library, params);
		
		ArrayList<XYPoint> points = getMatchingPoints(referenceRTInSecs, wideRTInSecs);
		
		RetentionTimeAlignmentInterface alignment=RetentionTimeFilter.getFilter(points, "reference", library.getName(), TwoDimensionalKDE.HIGHER_RESOLUTION);
		
		final List<RetentionTimeAlignmentInterface.AlignmentDataPoint> alignmentResults=alignment.plot(points, Optional.ofNullable(saveFileSeed));
		return new Pair<RetentionTimeAlignmentInterface, List<AlignmentDataPoint>>(alignment, alignmentResults);
	}

	public static ArrayList<XYPoint> getMatchingPoints(TObjectFloatHashMap<String> reference,
			TObjectFloatHashMap<String> wide) {
		ArrayList<XYPoint> points=new ArrayList<XYPoint>();
		
		wide.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String a, float b) {
				float alt=reference.get(a);
				if (reference.getNoEntryValue()!=alt) {
					// narrow first, then wide
					points.add(new XYPoint(alt, b));
				}
				return true;
			}
		});
		return points;
	}

	private static TObjectFloatHashMap<String> getRTs(LibraryInterface library, SearchParameters params) throws IOException, SQLException, DataFormatException {
		ArrayList<LibraryEntry> entries=library.getAllEntries(false, params.getAAConstants());
		TObjectFloatHashMap<String> rtInSec=new TObjectFloatHashMap<>();
		for (LibraryEntry entry : entries) {
			rtInSec.put(entry.getPeptideModSeq(), entry.getRetentionTime());
		}
		return rtInSec;
	}
}
