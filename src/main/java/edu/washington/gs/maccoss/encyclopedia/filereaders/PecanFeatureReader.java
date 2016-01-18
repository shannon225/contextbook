package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.AbstractPecanFragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringResultsToTSVConsumer;
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

public class PecanFeatureReader {
	public static ArrayList<LibraryEntry> parsePecanFeatures(File f, ArrayList<ScoredObject<String>> scoredPSMIDs, StripeFileInterface stripeFile, PecanScoringFactory factory) {
		final TObjectFloatHashMap<String> savedIDs=new TObjectFloatHashMap<String>();
		for (ScoredObject<String> psm : scoredPSMIDs) {
			savedIDs.put(psm.y, psm.x);
		}

		int cores=factory.getParameters().getNumberOfThreadsUsed();
		
		BlockingQueue<Map<String, String>> blockingQueue=new LinkedBlockingQueue<Map<String, String>>();
		TableParserProducer producer=new TableParserProducer(blockingQueue, f, "\t", cores);
		
		ConcurrentLinkedQueue<LibraryEntry> savedEntries=new ConcurrentLinkedQueue<LibraryEntry>();

		Thread producerThread=new Thread(producer);
		producerThread.start();

		Thread[] consumers=new Thread[cores];
		for (int i=0; i<consumers.length; i++) {
			PecanFeatureMuscle muscle=new PecanFeatureMuscle(savedIDs, stripeFile, factory, savedEntries);
			TableParserConsumer consumer=new TableParserConsumer(blockingQueue, muscle);
			consumers[i]=new Thread(consumer);
			consumers[i].start();
		}

		while (isAlive(consumers)) {
			String message=savedEntries.size()+" of "+scoredPSMIDs.size()+" peptides re-extracted";
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

	public static class PecanFeatureMuscle implements TableParserMuscle {
		private final TObjectFloatHashMap<String> savedIDs;
		private final StripeFileInterface stripeFile;

		private final PSMScorer scorer;
		private final SearchParameters params;
		private final PecanScoringFactory factory;

		private final ConcurrentLinkedQueue<LibraryEntry> savedEntries;

		public PecanFeatureMuscle(TObjectFloatHashMap<String> savedIDs, StripeFileInterface stripeFile, PecanScoringFactory factory, ConcurrentLinkedQueue<LibraryEntry> savedEntries) {
			this.savedIDs=new TObjectFloatHashMap<String>(savedIDs); // to guarantee immutability
			this.stripeFile=stripeFile;

			scorer=factory.getPecanScorer();
			params=factory.getParameters();
			this.factory=factory;
			this.savedEntries=savedEntries;
		}

		@Override
		public void processRow(Map<String, String> row) {
			String psmID=row.get("id");
			if (savedIDs.contains(psmID)) {
				int scanID=Integer.parseInt(row.get("ScanNr"));
				double precursorMZ=Double.parseDouble(row.get("precursorMz"));
				byte precursorCharge=PecanScoringResultsToTSVConsumer.getCharge(psmID);
				String peptideModSeq=PecanScoringResultsToTSVConsumer.getPeptideSequence(psmID);
				int copies=1;
				float retentionTime=Float.parseFloat(row.get("midTime"));
				float score=savedIDs.get(psmID);

				float duration=Float.parseFloat(row.get("duration"));

				try {
					ArrayList<Stripe> stripes=stripeFile.getStripes(precursorMZ, retentionTime-duration, retentionTime+duration, false);

					float bestDelta=Float.MAX_VALUE;
					Stripe bestStripe=null;
					for (Stripe stripe : stripes) {
						float delta=Math.abs(stripe.getScanStartTime()-retentionTime);
						if (delta<bestDelta) {
							bestDelta=delta;
							bestStripe=stripe;
						}
					}

					AbstractPecanFragmentationModel model=factory.getFragmentationModel(new FastaEntry(peptideModSeq), params.getAAConstants());
					LibraryEntry unitEntry=model.getUnitSpectrum(precursorCharge, params);
					PeakScores[] individualPeakScores=scorer.getIndividualPeakScores(unitEntry, bestStripe, true);

					TDoubleArrayList mzs=new TDoubleArrayList();
					TFloatArrayList intens=new TFloatArrayList();
					for (PeakScores scores : individualPeakScores) {
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
