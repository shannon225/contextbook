package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;
import gnu.trove.map.hash.TCharDoubleHashMap;
import junit.framework.TestCase;

public class PercolatorExecutorTest extends TestCase {
	public static void main(String[] args) throws Exception {
		File fastaFile=new File("/Users/searleb/Documents/vaneyk/control/Synthetic_Peptides.fasta");
		File featureFile=new File("/Users/searleb/Documents/vaneyk/control/only_synthesized_peptides_features.txt");
		File outputFile=new File("/Users/searleb/Documents/vaneyk/control/only_synthesized_peptides_concatenated_results.txt");
		File decoyFile=new File("/Users/searleb/Documents/vaneyk/control/only_synthesized_peptides_concatenated_decoy.txt");
		File outputProteinFile=new File("/Users/searleb/Documents/vaneyk/control/only_synthesized_peptides_concatenated_protein_results.txt");
		File decoyProteinFile=new File("/Users/searleb/Documents/vaneyk/control/only_synthesized_peptides_concatenated_protein_decoy.txt");
		PercolatorExecutionData percolatorFiles=new PercolatorExecutionData(featureFile, fastaFile, outputFile, decoyFile, outputProteinFile, decoyProteinFile, SearchParameterParser.getDefaultParametersObject());
		PercolatorExecutor e=new PercolatorExecutor(getDefaultPercolaterVersion(), percolatorFiles);
		BlockingQueue<OutputMessage> result=e.start();
             
		int outputlines=0;

		while (!e.isFinished()||!result.isEmpty()) {
			if (!result.isEmpty()) {
				OutputMessage data=result.take();
				if (data.isStdOutput()) {
					outputlines++;
				} else {
					System.out.println(data.getMessage());
				}
			} else {
				Thread.sleep(10);
			}
		}

		final AminoAcidConstants aaConstants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());

		System.out.println("total processed: "+outputlines);
		ArrayList<PercolatorPeptide> peptides=PercolatorReader.getPassingPeptidesFromTSV(outputFile, 0.01f, aaConstants, false).x;
		System.out.println("Peptides: "+peptides.size());
		ArrayList<PercolatorPeptide> decoys=PercolatorReader.getPassingPeptidesFromTSV(decoyFile, 0.01f, aaConstants, true).x;
		System.out.println("Decoys: "+decoys.size());
		
	}

	/**
	 * Used only in {@link #main}.
	 */
	static byte getDefaultPercolaterVersion() {
		switch (OSDetector.getOS()) {
			case WINDOWS:
				//TODO: issue #23: Percolator v3 fails silently with exit code 255 on some Windows machines
				return 2;
			default:
				return PercolatorExecutor.DEFAULT_VERSION_NUMBER;
		}
	}

	public void testParsePeptideSequence() {
		String peptideString="-.FNNFINDSLLEGAIDALKR.-";
		String parsed=PercolatorExecutor.parsePeptideSequence(peptideString);
		assertEquals("FNNFINDSLLEGAIDALKR", parsed);
	}

	public void testGetPercolatorVersionFromConsole() {
		String line = "Percolator version 3.01, Build Date May 23 2017 12:14:41";
		assertEquals("3.01", PercolatorExecutor.getPercolatorVersionFromOutput(line).orElse(null));

		line = "Percolator version 3.14.15, Build Date May 23 2017 12:14:41";
		assertEquals("3.14.15", PercolatorExecutor.getPercolatorVersionFromOutput(line).orElse(null));

		// Taken directly from executing src/main/resources/bin/percolator-v2-10.lin
		line = "Percolator version 2.09, Build Date Apr 15 2016 15:42:56";
		assertEquals("2.09", PercolatorExecutor.getPercolatorVersionFromOutput(line).orElse(null));

		line = "Percolator version 2, Build Date May 23 2017 12:14:41";
		assertEquals("2", PercolatorExecutor.getPercolatorVersionFromOutput(line).orElse(null));

		line = "Percolator version , Build Date May 23 2017 12:14:41";
		assertFalse(PercolatorExecutor.getPercolatorVersionFromOutput(line).isPresent());

		line = "Percolator version, Build Date May 23 2017 12:14:41";
		assertFalse(PercolatorExecutor.getPercolatorVersionFromOutput(line).isPresent());
	}

	public void testGetErrorMessageFromConsole() {
		assertEquals("bad allocation", PercolatorExecutor.getErrorMessage(new OutputMessage("Exception caught: bad allocation", false)));
		assertEquals("<error string>", PercolatorExecutor.getErrorMessage(new OutputMessage("Error : <error string>", false)));

		String msg = "anything that mentions a bad allocation";
		assertEquals(msg, PercolatorExecutor.getErrorMessage(new OutputMessage(msg, false)));
	}

	public static PercolatorExecutionData getPercolatorFiles(File featureFile, File fastaFile, SearchParameters parameters) throws IOException {
		File outputFile=File.createTempFile("percolator", ".txt");
		outputFile.deleteOnExit();
		File decoyFile=File.createTempFile("percolator", ".decoy.txt");
		decoyFile.deleteOnExit();
		File outputProteinFile=File.createTempFile("percolator", "protein.txt");
		outputProteinFile.deleteOnExit();
		File decoyProteinFile=File.createTempFile("percolator", ".protein_decoy.txt");
		decoyProteinFile.deleteOnExit();
		PercolatorExecutionData percolatorFiles=new PercolatorExecutionData(featureFile, fastaFile, outputFile, decoyFile, outputProteinFile, decoyProteinFile, parameters);
		return percolatorFiles;
	}
}
