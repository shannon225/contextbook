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

import edu.washington.gs.maccoss.encyclopedia.algorithms.DotProduct;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringResultsToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class SearchFeatureReader {
	public static ArrayList<LibraryEntry> parseSearchFeatures(File f, ArrayList<ScoredObject<String>> globalPassingPSMIDs, ArrayList<ScoredObject<String>> localPassingPSMIDs, StripeFileInterface stripeFile, SearchParameters parameters) {
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
		
		ConcurrentLinkedQueue<LibraryEntry> savedEntries=new ConcurrentLinkedQueue<LibraryEntry>();

		Thread producerThread=new Thread(producer);
		producerThread.start();

		Thread[] consumers=new Thread[cores];
		for (int i=0; i<consumers.length; i++) {
			SearchFeatureMuscle muscle=new SearchFeatureMuscle(savedIDs, stripeFile, parameters, savedEntries);
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
		
		ArrayList<LibraryEntry> entryList=new ArrayList<LibraryEntry>();
		for (LibraryEntry entry : savedEntries) {
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
		private final TObjectFloatHashMap<String> savedIDs;
		private final StripeFileInterface stripeFile;

		private final PSMScorer scorer;
		private final SearchParameters params;

		private final ConcurrentLinkedQueue<LibraryEntry> savedEntries;

		public SearchFeatureMuscle(TObjectFloatHashMap<String> savedIDs, StripeFileInterface stripeFile, SearchParameters parameters, ConcurrentLinkedQueue<LibraryEntry> savedEntries) {
			this.savedIDs=new TObjectFloatHashMap<String>(savedIDs); // to guarantee immutability
			this.stripeFile=stripeFile;

			scorer=new DotProduct(parameters.getFragmentTolerance());
			params=parameters;
			this.savedEntries=savedEntries;
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

				try {
					ArrayList<Stripe> stripes=stripeFile.getStripes(precursorMZ, retentionTime-duration, retentionTime+duration, false);

					FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
					LibraryEntry unitEntry=model.getUnitSpectrum(precursorCharge, params);

					float bestDelta=Float.MAX_VALUE;
					PeakScores[] bestScores=null;
					ArrayList<PeakScores[]> scoreList=new ArrayList<PeakScores[]>();
					for (Stripe stripe : stripes) {
						float delta=Math.abs(stripe.getScanStartTime()-retentionTime);
						PeakScores[] individualPeakScores=scorer.getIndividualPeakScores(unitEntry, stripe, true);
						scoreList.add(individualPeakScores);
						if (delta<bestDelta) {
							bestDelta=delta;
							bestScores=individualPeakScores;
						}
						if (peptideModSeq.equals("AAPQS[+80.0]PSVPK")) {
							System.out.println(peptideModSeq+": "+stripe.getSpectrumIndex()+"\t"+stripe.getScanStartTime()+"\t (target: "+retentionTime+")");
						}
					}
					
					int group=0;
					if (peptideModSeq.equals("AAPQS[+80.0]PSVPK")) {
						//System.out.println("ArrayList<PeakScores[]> scoreList=new ArrayList<PeakScores[]>();");
						for (PeakScores[] peakScores : scoreList) {
							//System.out.println("PeakScores[] peakScores=new PeakScores["+peakScores.length+"];");
							for (int i=0; i<peakScores.length; i++) {
								if (peakScores[i]!=null) {
									//System.out.println("peakScores["+i+"]=new PeakScores("+peakScores[i].getScore()+", "+peakScores[i].getTargetMass()+", "+peakScores[i].getDeltaMass()+");");
									System.out.println(group+"\t"+i+"\t"+peakScores[i].getScore()+"\t"+peakScores[i].getTargetMass()+"\t"+peakScores[i].getDeltaMass());
								}
							}
							//System.out.println("scoreList.add(peakScores);");
							group++;
						}
					}

					TDoubleArrayList mzs=new TDoubleArrayList();
					TFloatArrayList intens=new TFloatArrayList();
					for (PeakScores scores : bestScores) {
						if (scores!=null) {
							float peakScore=scores.getScore();
							if (peakScore>0) {
								mzs.add(scores.getTargetMass());
								intens.add(peakScore);
							}
						}
					}

					double[] massArray=mzs.toArray();
					float[] intensityArray=intens.toArray();
					LibraryEntry entry=new LibraryEntry(scanID, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray);
					savedEntries.add(entry);

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
		}
	};

}
