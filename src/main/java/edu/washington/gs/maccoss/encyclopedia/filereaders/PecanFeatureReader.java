package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
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
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
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

		PecanFeatureMuscle muscle=new PecanFeatureMuscle(savedIDs, stripeFile, factory);

		BlockingQueue<Map<String, String>> blockingQueue=new LinkedBlockingQueue<Map<String, String>>();
		TableParserProducer producer=new TableParserProducer(blockingQueue, f, "\t");
		TableParserConsumer consumer=new TableParserConsumer(blockingQueue, muscle);

		Thread producerThread=new Thread(producer);
		Thread consumerThread=new Thread(consumer);
		producerThread.start();
		consumerThread.start();

		try {
			producerThread.join();
			consumerThread.join();
		} catch (InterruptedException ie) {
			Logger.errorLine("Percolator reading interrupted!");
			Logger.errorException(ie);
		}
		
		return muscle.getSavedEntries();
	}

	public static class PecanFeatureMuscle implements TableParserMuscle {
		private final TObjectFloatHashMap<String> savedIDs;
		private final StripeFileInterface stripeFile;

		private final PSMScorer scorer;
		private final SearchParameters params;
		private final PecanScoringFactory factory;

		private final ArrayList<LibraryEntry> savedEntries=new ArrayList<LibraryEntry>();

		public PecanFeatureMuscle(TObjectFloatHashMap<String> savedIDs, StripeFileInterface stripeFile, PecanScoringFactory factory) {
			this.savedIDs=savedIDs;
			this.stripeFile=stripeFile;

			scorer=factory.getPecanScorer();
			params=factory.getParameters();
			this.factory=factory;
		}
		
		public ArrayList<LibraryEntry> getSavedEntries() {
			return savedEntries;
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

					AbstractPecanFragmentationModel model=factory.getFragmentationModel(peptideModSeq, params.getAAConstants());
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
					throw new EncyclopediaException("Error parsing Stripe file", ioe);
				} catch (SQLException sqle) {
					throw new EncyclopediaException("Error parsing Stripe file", sqle);
				} catch (DataFormatException dfe) {
					throw new EncyclopediaException("Error parsing Stripe file", dfe);
				}

			}
		}
	};

}
