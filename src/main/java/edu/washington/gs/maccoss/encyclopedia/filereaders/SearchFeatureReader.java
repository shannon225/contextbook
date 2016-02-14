package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.zip.DataFormatException;

import com.google.common.base.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.DotProduct;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringResultsToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoPermuter;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.set.hash.TFloatHashSet;

public class SearchFeatureReader {
	public static ArrayList<IntegratedLibraryEntry> parseSearchFeatures(File f, ArrayList<ScoredObject<String>> globalPassingPSMIDs, ArrayList<ScoredObject<String>> localPassingPSMIDs, StripeFileInterface stripeFile, Optional<LibraryFile> libraryFile, SearchParameters parameters) {
		HashSet<String> passingPeptideSequences=new HashSet<String>();
		for (ScoredObject<String> psm : globalPassingPSMIDs) {
			String peptideModSeq=PecanScoringResultsToTSVConsumer.getPeptideSequence(psm.y);
			passingPeptideSequences.add(peptideModSeq);
		}
		
		final TObjectFloatHashMap<String> savedIDs=new TObjectFloatHashMap<String>();
		for (ScoredObject<String> psm : localPassingPSMIDs) {
			String peptideModSeq=PecanScoringResultsToTSVConsumer.getPeptideSequence(psm.y);
			if (passingPeptideSequences.contains(peptideModSeq)) {
				savedIDs.put(psm.y, psm.x);
			}
		}

		int cores=parameters.getNumberOfThreadsUsed();
		
		BlockingQueue<Map<String, String>> blockingQueue=new LinkedBlockingQueue<Map<String, String>>();
		TableParserProducer producer=new TableParserProducer(blockingQueue, f, "\t", cores);
		
		ConcurrentLinkedQueue<IntegratedLibraryEntry> savedEntries=new ConcurrentLinkedQueue<IntegratedLibraryEntry>();

		Thread producerThread=new Thread(producer);
		producerThread.start();

		Thread[] consumers=new Thread[cores];
		for (int i=0; i<consumers.length; i++) {
			SearchFeatureMuscle muscle=new SearchFeatureMuscle(libraryFile, savedIDs, stripeFile, parameters, savedEntries);
			TableParserConsumer consumer=new TableParserConsumer(blockingQueue, muscle);
			consumers[i]=new Thread(consumer);
			consumers[i].start();
		}

		while (isAlive(consumers)) {
			String message=savedEntries.size()+" of "+globalPassingPSMIDs.size()+" peptides re-extracted";
			Logger.logLine(message+"...");
			try {
				Thread.sleep(500);
			} catch (InterruptedException ie) {
				ie.printStackTrace();
			}
		}

		try {
			producerThread.join();
			for (int i=0; i<consumers.length; i++) {
				consumers[i].join();
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("Percolator reading interrupted!");
			Logger.errorException(ie);
		}
		
		ArrayList<IntegratedLibraryEntry> entryList=new ArrayList<IntegratedLibraryEntry>();
		for (IntegratedLibraryEntry entry : savedEntries) {
			entryList.add(entry);
		}

		return entryList;
	}
	
	private static boolean isAlive(Thread[] consumers) {
		for (int i=0; i<consumers.length; i++) {
			if (consumers[i].isAlive()) return true;
		}
		return false;
	}

	public static class SearchFeatureMuscle implements TableParserMuscle {
		private final Optional<LibraryFile> library;
		private final TObjectFloatHashMap<String> savedIDs;
		private final StripeFileInterface stripeFile;
		private final boolean limitToQuantifiable;

		private final PSMScorer scorer;
		private final SearchParameters params;

		private final ConcurrentLinkedQueue<IntegratedLibraryEntry> savedEntries;

		public SearchFeatureMuscle(Optional<LibraryFile> library, TObjectFloatHashMap<String> savedIDs, StripeFileInterface stripeFile, SearchParameters parameters, ConcurrentLinkedQueue<IntegratedLibraryEntry> savedEntries) {
			this.library=library;
			this.savedIDs=new TObjectFloatHashMap<String>(savedIDs); // to guarantee immutability
			this.stripeFile=stripeFile;

			scorer=new DotProduct(parameters.getFragmentTolerance());
			params=parameters;
			this.savedEntries=savedEntries;
			
			this.limitToQuantifiable=true; //library.isPresent();
		}

		@Override
		public void processRow(Map<String, String> row) {
			String psmID=row.get("id");
			if (savedIDs.contains(psmID)) {
				int scanID=Integer.parseInt(row.get("ScanNr"));
				double precursorMZ=Double.parseDouble(row.get("precursorMz"));
				// FIXME need to get peptide charge from window
				byte precursorCharge=PecanScoringResultsToTSVConsumer.getCharge(psmID);
				String peptideModSeq=PecanScoringResultsToTSVConsumer.getPeptideSequence(psmID);
				int copies=1;
				String rtString=row.get("midTime");
				if (rtString==null) rtString=row.get("RTinMin");
				float retentionTime=Float.parseFloat(rtString)*60f;
				float score=savedIDs.get(psmID);

				String samplingTimeString=row.get("sampledTimes");
				float duration=samplingTimeString==null?(params.getExpectedPeakWidth()):Float.parseFloat(samplingTimeString);
				
				Optional<Triplet<double[], float[], Range>> spectrum=extractSpectrum(library, precursorMZ, precursorCharge, peptideModSeq, retentionTime, duration, limitToQuantifiable);
				if (params.isRunPhosphoLocalization()) {
					ArrayList<String> permutations=PhosphoPermuter.getPermutations(peptideModSeq, params.getAAConstants());
					if (permutations.size()==1) {
						System.out.println("single\t"+peptideModSeq);
					} else {
						boolean multiple=extractPhosphoForms(precursorMZ, precursorCharge, permutations, retentionTime);
						System.out.println("multiple\t"+peptideModSeq+"\t"+multiple);
					}
				}
				if (spectrum.isPresent()) {
					// FIXME need to not add duplicates!!!! for now just run SQL:
					// delete from entries where RowId not in (SELECT MIN(RowId) FROM entries GROUP BY PeptideModSeq, PrecursorCharge)
					Triplet<double[], float[], Range> spec=spectrum.get();
					IntegratedLibraryEntry entry=new IntegratedLibraryEntry(scanID, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, spec.x, spec.y, spec.z);
					if (limitToQuantifiable) {
						if (entry.getIonCount()<4||entry.getTIC()<1.0f) return; // FIXME HACKS GALORE
					}
					savedEntries.add(entry);
				}
			}
		}

		private boolean extractPhosphoForms(double precursorMZ, byte precursorCharge, ArrayList<String> peptideModSeqs, float retentionTime) {
			float duration=6*60f; // search for 6 minutes
			
			EncyclopediaOneScorer encyclopediaScorer=new EncyclopediaOneScorer(params, null); // not using aux scoring

			try {
				ArrayList<Stripe> stripes=stripeFile.getStripes(precursorMZ, retentionTime-duration, retentionTime+duration, false);
				
				TFloatArrayList allBestTimes=new TFloatArrayList();

				TFloatHashSet list=new TFloatHashSet();
				for (String peptideModSeq : peptideModSeqs) {
					FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
					LibraryEntry unitEntry=model.getUnitSpectrum(precursorCharge, params);	

					TFloatArrayList bestTimes=new TFloatArrayList();
					float bestScore=0.0f;
					TFloatFloatHashMap rtScoreMap=new TFloatFloatHashMap();
					for (Stripe spectrum : stripes) {
						//float score=encyclopediaScorer.score(unitEntry, spectrum);
						PeakScores[] individualPeakScores=encyclopediaScorer.getIndividualPeakScores(unitEntry, spectrum, false);
						float score=0.0f;
						for (int i=0; i<individualPeakScores.length; i++) {
							if (individualPeakScores[i]!=null) {
								score++;
							}
						}
						rtScoreMap.put(spectrum.getScanStartTime(), score);
						if (score>9.5f) { // at least 5 transitions
							if (bestScore<score) {
								bestScore=score;
								bestTimes.clear();
								bestTimes.add(spectrum.getScanStartTime());
							} else if (bestScore==score) {
								bestTimes.add(spectrum.getScanStartTime());
							}
						}
					}
					
					allBestTimes.addAll(bestTimes);
					
					EValueCalculator calculator=new EValueCalculator(rtScoreMap);
					//System.out.println(peptideModSeq+"\t"+calculator.getMaxRT()+"\t"+calculator.getNegLog10EValue()+"\t"+calculator.getMaxRawScore()); //FIXME
					list.add(calculator.getMaxRT());
				}
				
				if (allBestTimes.size()>1) {
					float range=allBestTimes.max()-allBestTimes.min();
					if (range>=60) {
						//System.out.println("multiple\t"+peptideModSeqs.get(0));
						return true;
					}
				} else if (allBestTimes.size()==1) {
					//System.out.println("single\t"+peptideModSeqs.get(0));
				}
				return false;

			} catch (IOException ioe) {
				Logger.errorLine("Error processing "+stripeFile.getFile().getName());
				throw new EncyclopediaException("Error parsing Stripe file", ioe);
			} catch (SQLException sqle) {
				Logger.errorLine("Error processing "+stripeFile.getFile().getName());
				throw new EncyclopediaException("Error parsing Stripe file", sqle);
			} catch (DataFormatException dfe) {
				Logger.errorLine("Error processing "+stripeFile.getFile().getName());
				throw new EncyclopediaException("Error parsing Stripe file", dfe);
			}
		}

		private Optional<Triplet<double[], float[], Range>> extractSpectrum(Optional<LibraryFile> library, double precursorMZ, byte precursorCharge, String peptideModSeq, float retentionTime, float duration, boolean limitToQuantifiable) {
			LibraryEntry unitEntry=null;
			if (library.isPresent()) {
				try {
					ArrayList<LibraryEntry> entries=library.get().getEntries(peptideModSeq, precursorCharge, false);
					if (entries.size()>0) {
						unitEntry=entries.get(0).toUnitSpectrum();
					} else {
						 // if library is ok but spectrum is not in library, just return null (don't quantify)
						return Optional.absent();
					}
					if (unitEntry.getIonCount()<4) { // FIXME HACKS GALORE
						 // if unit spectrum has fewer than 4 peaks, just return null (don't quantify)
						return Optional.absent();
					}
				} catch (IOException ioe) {
					Logger.errorLine("Error processing "+stripeFile.getFile().getName());
					throw new EncyclopediaException("Error parsing Stripe file", ioe);
				} catch (SQLException sqle) {
					Logger.errorLine("Error processing "+stripeFile.getFile().getName());
					throw new EncyclopediaException("Error parsing Stripe file", sqle);
				} catch (DataFormatException dfe) {
					Logger.errorLine("Error processing "+stripeFile.getFile().getName());
					throw new EncyclopediaException("Error parsing Stripe file", dfe);
				}
			}
			
			if (unitEntry==null) {
				FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
				unitEntry=model.getUnitSpectrum(precursorCharge, params);
			}
			
			return Optional.fromNullable(extractSpectrum(unitEntry, precursorMZ, peptideModSeq, retentionTime, duration, limitToQuantifiable));
		}

		public Triplet<double[], float[], Range> extractSpectrum(LibraryEntry unitEntry, double precursorMZ, String peptideModSeq, float retentionTime, float duration, boolean limitToQuantifiable) {
			try {
				ArrayList<Stripe> stripes=stripeFile.getStripes(precursorMZ, retentionTime-duration, retentionTime+duration, false);


				return quantifyPeptide(scorer, unitEntry, peptideModSeq, retentionTime, limitToQuantifiable, stripes);

			} catch (IOException ioe) {
				Logger.errorLine("Error processing "+stripeFile.getFile().getName());
				throw new EncyclopediaException("Error parsing Stripe file", ioe);
			} catch (SQLException sqle) {
				Logger.errorLine("Error processing "+stripeFile.getFile().getName());
				throw new EncyclopediaException("Error parsing Stripe file", sqle);
			} catch (DataFormatException dfe) {
				Logger.errorLine("Error processing "+stripeFile.getFile().getName());
				throw new EncyclopediaException("Error parsing Stripe file", dfe);
			}
		}

	};

	public static Triplet<double[], float[], Range> quantifyPeptide(PSMScorer scorer, LibraryEntry unitEntry, String peptideModSeq, float retentionTime, boolean limitToQuantifiable, ArrayList<Stripe> stripes) {
		float bestDelta=Float.MAX_VALUE;
		PeakScores[] bestScores=null;
		ArrayList<PeakScores[]> scoreList=new ArrayList<PeakScores[]>();
		TFloatArrayList retentionTimes=new TFloatArrayList();
		for (Stripe stripe : stripes) {
			retentionTimes.add(stripe.getScanStartTime());
			float delta=Math.abs(stripe.getScanStartTime()-retentionTime);
			PeakScores[] individualPeakScores=scorer.getIndividualPeakScores(unitEntry, stripe, true);
			scoreList.add(individualPeakScores);
			if (delta<bestDelta) {
				bestDelta=delta;
				bestScores=individualPeakScores;
			}
		}
		
		// no signal of any kind at retention time!
		if (bestScores==null) return null;
		
		TFloatArrayList[] traces=new TFloatArrayList[bestScores.length];
		for (int i=0; i<traces.length; i++) {
			traces[i]=new TFloatArrayList();
		}
		for (PeakScores[] peakScores : scoreList) {
			for (int i=0; i<peakScores.length; i++) {
				if (peakScores[i]!=null) {
					traces[i].add(peakScores[i].getScore());
				} else {
					traces[i].add(0.0f);
				}
			}
		}
		
		ArrayList<PeakScores> keptPeaks=new ArrayList<PeakScores>();
		ArrayList<float[]> chromatograms=new ArrayList<float[]>();
		for (int i=0; i<bestScores.length; i++) {
			if (bestScores[i]!=null&&bestScores[i].getScore()>0) {
				float[] chromatogram=traces[i].toArray();
				chromatogram=SkylineSGFilter.paddedSavitzkyGolaySmooth(chromatogram);
				chromatograms.add(chromatogram);
				keptPeaks.add(bestScores[i]);
			}
		}
		Triplet<float[], float[], Range> trio=TransitionRefiner.identifyTransitions(peptideModSeq, chromatograms, retentionTimes.toArray());
		float[] correlations=trio.x;
		float[] integrations=trio.y;
		
		TDoubleArrayList mzs=new TDoubleArrayList();
		TFloatArrayList intens=new TFloatArrayList();
		int count=0;
		int quantCount=0;
		float correlationThreshold=0.0f;//limitToQuantifiable?TransitionRefiner.quantitativeCorrelationThreshold:TransitionRefiner.identificationCorrelationThreshold;
		for (int i=0; i<keptPeaks.size(); i++) {
			PeakScores scores=keptPeaks.get(i);
			if (correlations[i]>=correlationThreshold) {
				quantCount++;
				float peakScore=scores.getScore();
				if (peakScore>0) {
					mzs.add(scores.getTargetMass());
					intens.add(integrations[i]);
					count++;
				}
			} else {
				mzs.add(scores.getTargetMass());
				intens.add(Float.MIN_VALUE);
			}
		}
		
		//System.out.println(peptideModSeq+"\t"+keptPeaks.size()+"\t"+count+"\t"+quantCount);

		double[] massArray=mzs.toArray();
		float[] intensityArray=intens.toArray();
		return new Triplet<double[], float[], Range>(massArray, intensityArray, trio.z);
	}
}
