package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

public class PeptideScoringResultsToTSVConsumer implements Runnable {
	private final BlockingQueue<PeptideScoringResult> resultsQueue;
	private final PrintWriter writer;

	public PeptideScoringResultsToTSVConsumer(File outputFile, BlockingQueue<PeptideScoringResult> resultsQueue) {
		this.resultsQueue=resultsQueue;
		try {
			writer=new PrintWriter(outputFile, "UTF-8");
		} catch (FileNotFoundException e) {
			throw new EncyclopediaException("Error setting up output file: "+outputFile.getAbsolutePath(), e);
		} catch (UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error setting up output file: "+outputFile.getAbsolutePath(), e);
		}
	}
	
	public void close() {
		writer.flush();
		writer.close();
	}
	
	@Override
	public void run() {
		boolean printedHeader=false; 
		try {
			while (true) {
				PeptideScoringResult result=resultsQueue.take();
				if (PeptideScoringResult.POISON_RESULT==result) break;
				if (!printedHeader) {
					writer.println("id\tTD\tScanNr\ttopx\tpeakBGScore\tdeltaCn\trawScore\tpeakSimilarity\tweightedRawScore\tnumAboveThresholdMatches\tnumMatches\taverageAbsPPM\taveragePPM\tisotopeDotProduct\tnumAboveThresholdPeakIons\tnumPeakIons\tmidTime\tpepLength\tcharge2\tcharge3\tsequence\tannotation");
					printedHeader=true;
				}
				
				LibraryEntry peptide=result.getEntry();
				int rank=1;
				
				float secondScore=0.0f;
				if (result.getGoodStripes().size()>1) {
					Pair<ScoredObject<Stripe>, float[]> second=result.getGoodStripes().get(1);
					secondScore=second.x.x;
				}				
				for (Pair<ScoredObject<Stripe>, float[]> goodStripe : result.getGoodStripes()) {
					float primaryScore=goodStripe.x.x;
					Stripe stripe=goodStripe.x.y;
					float[] auxScores=goodStripe.y;
					
					if (rank<=3) {
						float deltaCn=secondScore<=0?1.0f:(primaryScore-secondScore)/secondScore;
						writer.print(peptide.getPeptideModSeq()+"+"+peptide.getPrecursorCharge()+"\t"+(peptide.isDecoy()?-1:1)+"\t"+stripe.getSpectrumIndex()+"\t"+rank+"\t"+primaryScore+"\t"+deltaCn);
						for (float s : auxScores) {
							writer.print("\t"+s);
						}
						String sequence="-."+peptide.getPeptideSeq()+".-";
						
						String annotation=stripe.getSpectrumName();
						writer.print("\t"+stripe.getScanStartTime()+"\t"+peptide.getPeptideSeq().length()+"\t"+(peptide.getPrecursorCharge()==2?1:0)+"\t"+(peptide.getPrecursorCharge()==3?1:0)+"\t"+sequence+"\t"+annotation);
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
