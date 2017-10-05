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
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
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
		ArrayList<PercolatorPeptide> peptides=PercolatorReader.getPassingPeptidesFromTSV(outputFile, 0.01f, false).x;
		System.out.println("Peptides: "+peptides.size());
		ArrayList<PercolatorPeptide> decoys=PercolatorReader.getPassingPeptidesFromTSV(decoyFile, 0.01f, true).x;
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

		File outputFile=File.createTempFile("percolator", ".txt");
		outputFile.deleteOnExit();

		File decoyFile=File.createTempFile("percolator", ".decoy.txt");
		decoyFile.deleteOnExit();

		Files.copy(is, featureFile.toPath(), StandardCopyOption.REPLACE_EXISTING);

		Pair<ArrayList<PercolatorPeptide>, Float> pair=PercolatorExecutor.executePercolatorTSV(getDefaultPercolaterVersion(), featureFile, outputFile, decoyFile, 0.01f);
		assertEquals(405, pair.x.size());
		assertEquals(0.348315f, pair.y, 0.001f);
		
		pair=PercolatorReader.getPassingPeptidesFromTSV(outputFile, 0.01f, false);
		assertEquals(405, pair.x.size());
		assertEquals(0.348315f, pair.y, 0.001f);
		
		Pair<ArrayList<PercolatorPeptide>, Float> decoyPair=PercolatorReader.getPassingPeptidesFromTSV(decoyFile, 0.01f, true);
		assertEquals(3, decoyPair.x.size());
		assertEquals(0.0f, decoyPair.y, 0.001f);
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
