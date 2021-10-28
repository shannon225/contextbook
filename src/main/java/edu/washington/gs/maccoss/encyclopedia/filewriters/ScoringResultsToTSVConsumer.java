package edu.washington.gs.maccoss.encyclopedia.filewriters;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.RescoredPeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.RescoredSpectrumScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

import java.io.File;
import java.util.HashSet;
import java.util.concurrent.BlockingQueue;

public class ScoringResultsToTSVConsumer extends AbstractScoringResultsToTSVConsumer {
	private final String[] scoreNames;
	private final SearchParameters params;

	public ScoringResultsToTSVConsumer(File outputFile, StripeFileInterface diaFile, String[] scoreNames, BlockingQueue<AbstractScoringResult> resultsQueue, SearchParameters params) {
		super(outputFile, diaFile, resultsQueue);
		this.scoreNames = scoreNames;
		this.params=params;
	}

	@Override
	public void run() {
		boolean printedHeader=false; 
		try {
			while (true) {
				AbstractScoringResult result=resultsQueue.take();
				
				if (AbstractScoringResult.POISON_RESULT==result) break;
				if (!printedHeader) {
					writer.print("id\tTD\tScanNr\t");
					for (String name : scoreNames) {
						writer.print(name);
						writer.print('\t');
					}
					
					if (result instanceof RescoredPeptideScoringResult||result instanceof RescoredSpectrumScoringResult) {
						writer.print("deltaRT\t");//discriminantScore\t");
					}
					writer.print("numMissedCleavage\tpepLength\tcharge1\tcharge2\tcharge3\tcharge4\tprecursorMz\tprecursorMass\tRTinMin\tsequence\tprotein");
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

	protected void processResult(AbstractScoringResult result) {
		if (!result.hasScoredResults()) return;
		
		LibraryEntry peptide=result.getEntry();
		int rank=1;

		for (Pair<ScoredObject<FragmentScan>, float[]> goodStripe : result.getGoodMSMSCandidates()) {
			numberProcessed++;

			FragmentScan stripe=goodStripe.x.y;
			float[] auxScores=goodStripe.y;


			if (stripe!=null) {
				String psmID= PercolatorPeptide.getPSMID(peptide, stripe.getScanStartTime(), diaFile);

				writer.print(psmID);
				writer.print("\t"+(peptide.isDecoy()?-1:1));
				writer.print("\t"+stripe.getSpectrumIndex());

				for (int i=0; i<auxScores.length; i++) {
					writer.print('\t');
					writer.print(auxScores[i]);
				}

				writer.print("\t"+params.getEnzyme().getNumMissedCleavages(peptide.getPeptideSeq()));
				writer.print("\t"+peptide.getPeptideSeq().length());
				writer.print("\t"+(peptide.getPrecursorCharge()==1?1:0));
				writer.print("\t"+(peptide.getPrecursorCharge()==2?1:0));
				writer.print("\t"+(peptide.getPrecursorCharge()==3?1:0));
				writer.print("\t"+(peptide.getPrecursorCharge()>=4?1:0));
				writer.print("\t"+peptide.getPrecursorMZ());
				writer.print("\t"+(peptide.getPrecursorMZ()*peptide.getPrecursorCharge()-MassConstants.protonMass*peptide.getPrecursorCharge()));
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
			if (rank>2) break;
		}
	}

	public String[] getScoreNames() {
		return scoreNames;
	}
}
