package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;
import junit.framework.TestCase;

public class PercolatorExecutorTest extends TestCase {
	public static void main(String[] args) throws Exception {
		File featureFile=new File("/Users/searleb/Documents/school/projects/test/02may2016_yeast_deep_dia_01.mzML.encyclopedia.txt.features.txt");
		File outputFile=new File("/Users/searleb/Documents/school/projects/test/02may2016_yeast_deep_dia_01.mzML.encyclopedia.txt.output.txt");
		File decoyFile=new File("/Users/searleb/Documents/school/projects/test/02may2016_yeast_deep_dia_01.mzML.encyclopedia.txt.decoy.txt");
		PercolatorExecutor e=new PercolatorExecutor(featureFile, outputFile, decoyFile, getDefaultPercolaterVersion(), false);
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
		System.out.println("total processed: "+outputlines);
		ArrayList<PercolatorPeptide> peptides=PercolatorReader.getPassingPeptidesFromTSV(outputFile, 0.01f);
		System.out.println("Peptides: "+peptides.size());
		ArrayList<PercolatorPeptide> decoys=PercolatorReader.getPassingPeptidesFromTSV(decoyFile, 0.01f);
		System.out.println("Decoys: "+decoys.size());
		
	}

	public void testParsePeptideSequence() {
		String peptideString="-.FNNFINDSLLEGAIDALKR.-";
		String parsed=PercolatorExecutor.parsePeptideSequence(peptideString);
		assertEquals("FNNFINDSLLEGAIDALKR", parsed);
	}

	public void testPercolatorExecutor() throws Exception {
		InputStream is=getClass().getResourceAsStream("/pecan.feature.txt");
		File featureFile=File.createTempFile("pecan", ".feature");
		featureFile.deleteOnExit();

		File outputFile=File.createTempFile("percolator", ".xml");
		outputFile.deleteOnExit();

		File decoyFile=File.createTempFile("percolator", ".decoy.xml");
		decoyFile.deleteOnExit();

		Files.copy(is, featureFile.toPath(), StandardCopyOption.REPLACE_EXISTING);

		PercolatorExecutor e=new PercolatorExecutor(featureFile, outputFile, decoyFile, getDefaultPercolaterVersion(), true);
		BlockingQueue<OutputMessage> result=e.start();

		int outputlines=0;

		while (!e.isFinished()||!result.isEmpty()) {
			if (!result.isEmpty()) {
				OutputMessage data=result.take();
				if (data.isStdOutput()) {
					outputlines++;
//					Logger.logLine("[percolator:stdout]" + data.getMessage());
				} else {
					// ensure that any error messages are written to the console for debugging
					Logger.logLine("[percolator:stderr]" + data.getMessage());
				}
			} else {
				Thread.sleep(10);
			}
		}

		assertEquals("Non-zero exit code!", 0, e.getResultCode());

		assertEquals("Wrong number of spectra above 1% FDR!", 712, outputlines-1); // number of spectra above 1% FDR (-1 for header)
	}

	//TODO: issue #23: Percolator v3 fails silently with exit code 255 on some Windows machines
	private static byte getDefaultPercolaterVersion() {
		switch (OSDetector.getOS()) {
			case WINDOWS:
				return 2;
			default:
				return PercolatorExecutor.DEFAULT_VERSION_NUMBER;
		}
	}
}
