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
					writer.println("id\tTD\tScanNr\ttopx\tpeakBGScore\tdeltaCn\ttraceNumAboveThresholdIons\ttraceNumIons\tmidTime\t"
							+ "peakRawScore\tpeakSimilarity\tpeakWeightedRawScore\tpeakNumAboveThresholdMatches\tpeakNumMatches\tpeakAverageAbsPPM\tpeakAveragePPM\tpeakIsotopeDotProduct\t"
							+ "midRawScore\tmidSimilarity\tmidWeightedRawScore\tmidNumAboveThresholdIons\tmidNumIons\tmidAbsPPM\tmidPPM\tmidIsotopeDotProduct\t"
							+ "pepLength\tcharge2\tcharge3\tsequence\tannotation");
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
						writer.print((peptide.isDecoy()?"decoy":"")+peptide.getPeptideModSeq()+"+"+peptide.getPrecursorCharge()+"\t"+(peptide.isDecoy()?-1:1)+"\t"+stripe.getSpectrumIndex()+"\t"+rank+"\t"+primaryScore+"\t"+deltaCn);
						for (float s : auxScores) {
							writer.print("\t"+s);
						}
						String sequence="-."+peptide.getPeptideSeq()+".-";
						
						String annotation=stripe.getSpectrumName();
						writer.print("\t"+peptide.getPeptideSeq().length()+"\t"+(peptide.getPrecursorCharge()==2?1:0)+"\t"+(peptide.getPrecursorCharge()==3?1:0)+"\t"+sequence+"\t"+annotation);
						writer.println();
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
