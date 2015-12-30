package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

public class PecanScoringResultsToTSVConsumer implements PeptideScoringResultsConsumer {
	private final BlockingQueue<PeptideScoringResult> resultsQueue;
	private final PrintWriter writer;
	private volatile int numberProcessed=0;
	private final int numberOfPeaksPerPeptide;

	public PecanScoringResultsToTSVConsumer(File outputFile, BlockingQueue<PeptideScoringResult> resultsQueue, int numberOfPeaksPerPeptide) {
		this.resultsQueue=resultsQueue;
		this.numberOfPeaksPerPeptide=numberOfPeaksPerPeptide;
		try {
			writer=new PrintWriter(outputFile, "UTF-8");
		} catch (FileNotFoundException e) {
			throw new EncyclopediaException("Error setting up output file: "+outputFile.getAbsolutePath(), e);
		} catch (UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error setting up output file: "+outputFile.getAbsolutePath(), e);
		}
	}

	@Override
	public void close() {
		writer.flush();
		writer.close();
	}

	@Override
	public int getNumberProcessed() {
		return numberProcessed;
	}
	
	@Override
	public void run() {
		boolean printedHeader=false; 
		try {
			while (true) {
				PeptideScoringResult result=resultsQueue.take();
				if (PeptideScoringResult.POISON_RESULT==result) break;
				if (!printedHeader) {
					// Percolator always assumes linux line endings!
					writer.print("id\tTD\tScanNr\ttopx\tpeakBGScore\tdeltaCn\t"
							+ "peakAvgIdotp\tpeakMaxIdotp\tpeakScore\tpeakWScore\tmidIons\tpeakMassErrMean\tpeakMassErrVar\tprecursorMassErrMean\t"
							+ "precursorMassErrVar\tpeakSimilarity\tduration\tmidTime\t"
							+ "pepLength\tcharge2\tcharge3\tsequence\tannotation\n");
					printedHeader=true;
				}
				LibraryEntry peptide=result.getEntry();
				int rank=1;
				
				float firstScore=0.0f;
				float secondScore=0.0f;

				if (result.getGoodStripes().size()>0) {
					Pair<ScoredObject<Stripe>, float[]> first=result.getGoodStripes().get(0);
					firstScore=first.x.x;
				}
				if (result.getGoodStripes().size()>1) {
					Pair<ScoredObject<Stripe>, float[]> second=result.getGoodStripes().get(1);
					secondScore=second.x.x;
				}
				for (Pair<ScoredObject<Stripe>, float[]> goodStripe : result.getGoodStripes()) {
					numberProcessed++;
					
					float primaryScore=goodStripe.x.x;
					Stripe stripe=goodStripe.x.y;
					float[] auxScores=goodStripe.y;
					
					if (rank<=numberOfPeaksPerPeptide) {
						float deltaCn=firstScore<=0?0.0f:Math.min(1.0f, (primaryScore-secondScore)/firstScore); // if secondScore<0 then deltaCn can be >1, so protect against that
						writer.print((peptide.isDecoy()?"decoy":"")+peptide.getPeptideModSeq()+"+"+peptide.getPrecursorCharge()+"\t"+(peptide.isDecoy()?-1:1)
								+"\t"+stripe.getSpectrumIndex()+"\t"+rank+"\t"+auxScores[16]+"\t"+deltaCn);

						/*
						 * 0) traceNumAboveThresholdIons  
						 * 1) traceNumIons  
						 * 2) midTime  
						 * 3) peakRawScore  
						 * 4) peakSimilarity  
						 * 5) peakWeightedRawScore  
						 * 6) peakNumAboveThresholdMatches  
						 * 7) peakNumMatches  
						 * 8) peakAverageAbsPPM  
						 * 9) peakAveragePPM  
						 * 10) peakIsotopeDotProduct  
						 * 11) fragmentDeltaMassAverage
						 * 12) fragmentDeltaMassVariance
						 * 13) duration
						 * 14) tpeakMaxIdotp
						 * 15) varPPM  
						 * 16) bgSubScore
						 * 
						 * peakAvgIdotp	peakMaxIdotp	peakScore	peakWScore	peakIons	peakMassErrMean	peakMassErrVar	precursorMassErrMean	
						 * precursorMassErrVar	peakSimilarity	duration	midTime
						 */

						writer.print("\t"+auxScores[10]); //peakIsotopeDotProduct
						writer.print("\t"+auxScores[14]); //tpeakMaxIdotp
						writer.print("\t"+auxScores[3]);  //peakRawScore
						writer.print("\t"+auxScores[5]);  //peakWeightedRawScore
						writer.print("\t"+auxScores[0]);  //traceNumAboveThresholdIons
						writer.print("\t"+auxScores[11]); //fragmentDeltaMassAverage
						writer.print("\t"+auxScores[12]); //fragmentDeltaMassVariance
						writer.print("\t"+auxScores[9]);  //peakAveragePPM
						writer.print("\t"+auxScores[15]); //varPPM
						writer.print("\t"+auxScores[4]);  //peakSimilarity
						writer.print("\t"+auxScores[13]); //duration
						writer.print("\t"+auxScores[2]);  //midTime
						String sequence="-."+peptide.getPeptideSeq()+".-";

						String annotation=stripe.getSpectrumName();
						writer.print("\t"+peptide.getPeptideSeq().length()+"\t"+(peptide.getPrecursorCharge()==2?1:0)+"\t"+(peptide.getPrecursorCharge()==3?1:0)+"\t"+sequence+"\t"+annotation);

						// Percolator always assumes linux line endings!
						writer.print("\n");
					}
					rank++;
					if (rank>3) break;
				}
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("DIA writing interrupted!");
			Logger.errorException(ie);
		}
	}
}
