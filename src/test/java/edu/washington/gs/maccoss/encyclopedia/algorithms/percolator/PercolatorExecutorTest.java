package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;
import junit.framework.TestCase;

public class PercolatorExecutorTest extends TestCase {
	
	public static void main(String[] args) throws Exception {
		File featureFile=new File("/Users/searleb/Documents/school/projects/pecandata/DIA_1xGFP_20x20mz_500to900_rep1.mzML.pecan.txt.features.txt");
		PercolatorExecutor e=new PercolatorExecutor(featureFile);
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
		Files.copy(is, featureFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
		
		PercolatorExecutor e=new PercolatorExecutor(featureFile);
		BlockingQueue<OutputMessage> result=e.start();
		
		int outputlines=0;

		while (!e.isFinished()||!result.isEmpty()) {
			if (!result.isEmpty()) {
				OutputMessage data=result.take();
				if (data.isStdOutput()) {
					outputlines++;
				}
			} else {
				Thread.sleep(10);
			}
		}
		assertEquals(712, outputlines-1); // number of spectra above 1% FDR (-1 for header)
	}
}
