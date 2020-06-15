package edu.washington.gs.maccoss.encyclopedia.filewriters;

import com.github.davidmoten.bigsorter.Serializer;
import com.github.davidmoten.bigsorter.Sorter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.RescoredPeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.io.FileUtils;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.Comparator;
import java.util.HashSet;
import java.util.concurrent.BlockingQueue;

public class ScoringResultsToTSVConsumer implements PeptideScoringResultsConsumer {
	private final OS os=OSDetector.getOS();

	private final File outputFile, tmpFile;
	private final StripeFileInterface diaFile;
	private final String[] scoreNames;
	private final BlockingQueue<PeptideScoringResult> resultsQueue;
	private final int numberOfPeaksPerPeptide;

	private final PrintWriter writer;

	private volatile int numberProcessed=0;

	public ScoringResultsToTSVConsumer(File outputFile, StripeFileInterface diaFile, String[] scoreNames, BlockingQueue<PeptideScoringResult> resultsQueue, int numberOfPeaksPerPeptide) {
		this.outputFile = outputFile;
		this.tmpFile = new File(outputFile.getAbsolutePath() + ".unsorted");
		this.diaFile=diaFile;
		this.resultsQueue=resultsQueue;
		this.numberOfPeaksPerPeptide=numberOfPeaksPerPeptide;
		this.scoreNames=scoreNames;

		try {
			writer=new PrintWriter(this.tmpFile, "UTF-8");
			System.out.println("Constructing writer for "+ this.tmpFile.getAbsolutePath());
		} catch (FileNotFoundException e) {
			throw new EncyclopediaException("Error setting up output file: "+ this.tmpFile.getAbsolutePath(), e);
		} catch (UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error setting up output file: "+ this.tmpFile.getAbsolutePath(), e);
		}
	}

	public String[] getScoreNames() {
		return scoreNames;
	}

	@Override
	public BlockingQueue<PeptideScoringResult> getResultsQueue() {
		return resultsQueue;
	}

	@Override
	public void close() {
		writer.flush();
		writer.close();

		Logger.logLine("Sorting results into "+ this.outputFile.getAbsolutePath());

		try {
			final Serializer<CSVRecord> serializer = Serializer.csv(
					CSVFormat
							.DEFAULT
							.withDelimiter('\t')
							.withRecordSeparator(getLineSeparator())
							.withFirstRecordAsHeader(),
					StandardCharsets.UTF_8
			);

			// compare rows on column zero (ID) using lexicographic string ordering
			final Comparator<CSVRecord> comparator = Comparator.comparing(
					record -> record.get(0),
					Comparator.naturalOrder()
			);

			Sorter
					.serializer(serializer)
					.comparator(comparator)
					.input(this.tmpFile)
					.output(this.outputFile)
					.tempDirectory(this.outputFile.getParentFile()) // always sort in the target directory
					.maxItemsPerFile(100000) // 100k lines per file; this controls memory usage
					.sort();
		} catch (UncheckedIOException exception) {
			Logger.errorLine("Caught IO exception sorting TSV output; failing!");
			Logger.errorException(exception);
			throw exception;
		} finally {
			// unconditionally remove the unsorted temp file
			System.out.println("Removing temp file "+ this.tmpFile.getAbsolutePath());
			FileUtils.deleteQuietly(this.tmpFile);
		}
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

	private String getLineSeparator() {
		switch (os) {
			case MAC:
				return "\n";
			default:
				return System.lineSeparator();
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
}
