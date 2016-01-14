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
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

public class PecanScoringResultsToTSVConsumer implements PeptideScoringResultsConsumer {
	private final OS os=OSDetector.getOS();
	
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
					writer.print("id\tTD\tScanNr\ttopN\trank\tpeakZScore\tpeakCalibratedScore\tdeltaSn\t"
							+ "avgIdotp\tmidIdotp\tpeakScore\tpeakWeightedScore\tNCI\tCIMassErrMean\tCIMassErrVar\tprecursorMassErrMean\t"
							+ "precursorMassErrVar\tpeakSimilarity\tsampledTimes\tmidTime\tspectraNorm\t"
							+ "pepLength\tcharge2\tcharge3\tprecursorMz\tsequence\tannotation");
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
						String psmID = getPSMID(peptide);

						/*
						 * 0) traceNumAboveThresholdIons  
						 * 1) traceNumIons  
						 * 2) midTime  
						 * 3) peakRawScore  
						 * 4) peakSimilarity  
						 * 5) peakEuclideanDistance
						 * 6) peakWeightedRawScore  
						 * 7) peakNumAboveThresholdMatches  
						 * 8) peakNumMatches  
						 * 9) peakAverageAbsPPM  
						 * 10) peakAveragePPM  
						 * 11) peakIsotopeDotProduct  
						 * 12) fragmentDeltaMassAverage
						 * 13) fragmentDeltaMassVariance
						 * 14) duration
						 * 15) tpeakMaxIdotp
						 * 16) tpeakMidIdotp
						 * 17) varPPM  
						 * 18) bgSubScore
						 * 19) zScore
						 * 20) rank
						 * 
						 * peakAvgIdotp	peakMaxIdotp	peakScore	peakWScore	peakIons	peakMassErrMean	peakMassErrVar	precursorMassErrMean	
						 * precursorMassErrVar	peakSimilarity	duration	midTime
						 */

						writer.print(psmID);
						writer.print("\t"+(peptide.isDecoy()?-1:1));
						writer.print("\t"+stripe.getSpectrumIndex());
						writer.print("\t"+rank);
						writer.print("\t"+auxScores[20]);
						writer.print("\t"+auxScores[19]);
						writer.print("\t"+auxScores[18]);
						writer.print("\t"+deltaCn);
						writer.print("\t"+auxScores[11]); //peakIsotopeDotProduct
						writer.print("\t"+auxScores[16]); //tpeakMidIdotp
						writer.print("\t"+auxScores[3]);  //peakRawScore
						writer.print("\t"+auxScores[6]);  //peakWeightedRawScore
						writer.print("\t"+auxScores[0]);  //traceNumAboveThresholdIons
						writer.print("\t"+auxScores[12]); //fragmentDeltaMassAverage
						writer.print("\t"+auxScores[13]); //fragmentDeltaMassVariance
						writer.print("\t"+auxScores[10]);  //peakAveragePPM
						writer.print("\t"+auxScores[17]); //varPPM
						writer.print("\t"+auxScores[4]);  //peakSimilarity
						writer.print("\t"+auxScores[14]); //duration
						writer.print("\t"+auxScores[2]);  //midTime
						writer.print("\t"+auxScores[5]);  //peakEuclideanDistance

						writer.print("\t"+peptide.getPeptideSeq().length());
						writer.print("\t"+(peptide.getPrecursorCharge()==2?1:0));
						writer.print("\t"+(peptide.getPrecursorCharge()==3?1:0));
						writer.print("\t"+peptide.getPrecursorMZ());

						String sequence="-."+peptide.getPeptideSeq()+".-";
						String annotation;
						if (peptide instanceof PecanLibraryEntry) {
							annotation=((PecanLibraryEntry)peptide).getAccession();
						} else {
							annotation=stripe.getSpectrumName();
						}
						writer.print("\t"+sequence);
						writer.print("\t"+annotation);

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
		} catch (InterruptedException ie) {
			Logger.errorLine("DIA writing interrupted!");
			Logger.errorException(ie);
		}
	}

	public static String getPSMID(LibraryEntry peptide) {
		return (peptide.isDecoy()?"decoy":"")+peptide.getPeptideModSeq()+"+"+peptide.getPrecursorCharge();
	}

	public static boolean isPSMIDDecoy(String psmID) {
		return psmID.startsWith("decoy");
	}
	
	public static String getPeptideSequence(String psmID) {
		if (psmID.startsWith("decoy")) {
			psmID=psmID.substring(5);
		}
		return psmID.substring(0, psmID.indexOf('+'));
	}
	public static byte getCharge(String psmID) {
		return Byte.parseByte(psmID.substring(psmID.indexOf('+')+1));
	}
}
