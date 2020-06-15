package edu.washington.gs.maccoss.encyclopedia.filewriters;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.RescoredPeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

import java.io.File;
import java.util.HashSet;
import java.util.concurrent.BlockingQueue;

public class ScoringResultsToTSVConsumer extends AbstractScoringResultsToTSVConsumer {
	private final String[] scoreNames;
	private final int numberOfPeaksPerPeptide;

	public ScoringResultsToTSVConsumer(File outputFile, StripeFileInterface diaFile, String[] scoreNames, BlockingQueue<PeptideScoringResult> resultsQueue, int numberOfPeaksPerPeptide) {
		super(outputFile, diaFile, resultsQueue);
		this.scoreNames = scoreNames;
		this.numberOfPeaksPerPeptide=numberOfPeaksPerPeptide;
	}

	@Override
	public void run() {
		boolean printedHeader=false; 
		try {
			while (true) {
				PeptideScoringResult result=resultsQueue.take();
				
				if (PeptideScoringResult.POISON_RESULT==result) break;
				if (!printedHeader) {
					writer.print("id\tTD\tScanNr\ttopN\tdeltaCN\t");
					for (String name : scoreNames) {
						writer.print(name);
						writer.print('\t');
					}
					if (result instanceof RescoredPeptideScoringResult) {
						writer.print("deltaRT\t");//discriminantScore\t");
					}
					writer.print("pepLength\tcharge2\tcharge3\tprecursorMz\tRTinMin\tsequence\tprotein");
					// Percolator assumes linux line endings on Mac!
					switch (os) {
						case MAC:
							writer.print("\n");
							break;

						default:
							writer.println();
							break;
					}
					printedHeader=true;
				}
				processResult(result);
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("DIA writing interrupted!");
			Logger.errorException(ie);
		}
	}

	protected void processResult(PeptideScoringResult result) {
		LibraryEntry peptide=result.getEntry();
		int rank=1;

		float firstScore=0.0f;
		float secondScore=0.0f;

		if (result.getGoodStripes().size()>0) {
			Pair<ScoredObject<FragmentScan>, float[]> first=result.getGoodStripes().get(0);
			firstScore=first.x.x;
		}
		if (result.getGoodStripes().size()>1) {
			Pair<ScoredObject<FragmentScan>, float[]> second=result.getGoodStripes().get(1);
			secondScore=second.x.x;
		}

		for (Pair<ScoredObject<FragmentScan>, float[]> goodStripe : result.getGoodStripes()) {
			numberProcessed++;

			float primaryScore=goodStripe.x.x;
			FragmentScan stripe=goodStripe.x.y;
			float[] auxScores=goodStripe.y;


			if (stripe!=null&&rank<=numberOfPeaksPerPeptide) {
				float deltaCn;
				if (!Float.isNaN(firstScore)||!Float.isNaN(primaryScore)||!Float.isNaN(secondScore)||firstScore<=0) {
					deltaCn=0.0f;
				} else {
					deltaCn=Math.min(1.0f, (primaryScore-secondScore)/firstScore); // if secondScore<0 then deltaCn can be >1, so protect against that
				}
				String psmID= PercolatorPeptide.getPSMID(peptide, stripe.getScanStartTime(), diaFile);

				writer.print(psmID);
				writer.print("\t"+(peptide.isDecoy()?-1:1));
				writer.print("\t"+stripe.getSpectrumIndex());
				writer.print("\t"+rank);
				writer.print("\t"+deltaCn);

				for (int i=0; i<auxScores.length; i++) {
					writer.print('\t');
					writer.print(auxScores[i]);
				}

				writer.print("\t"+peptide.getPeptideSeq().length());
				writer.print("\t"+(peptide.getPrecursorCharge()==2?1:0));
				writer.print("\t"+(peptide.getPrecursorCharge()==3?1:0));
				writer.print("\t"+peptide.getPrecursorMZ());
				writer.print("\t"+stripe.getScanStartTime()/60f);

				String sequence="-."+peptide.getPeptideModSeq()+".-";
				writer.print("\t"+sequence);

				HashSet<String> accessions=peptide.getAccessions();
				writer.print("\t"+ PSMData.accessionsToString(accessions));

				// Percolator assumes linux line endings on Mac!
				switch (os) {
					case MAC:
						writer.print("\n");
						break;

					default:
						writer.println();
						break;
				}
			}
			rank++;
			if (rank>3) break;
		}
	}

	public String[] getScoreNames() {
		return scoreNames;
	}
}
