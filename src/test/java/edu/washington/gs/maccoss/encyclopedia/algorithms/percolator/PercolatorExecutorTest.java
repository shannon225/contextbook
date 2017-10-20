package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;
import junit.framework.TestCase;

public class PercolatorExecutorTest extends TestCase {
	public static void main(String[] args) throws Exception {
		File fastaFile=new File("/Users/searleb/Documents/projects/phosphopedia/sp_iso_HUMAN_4.9.2015_UP000005640_plusReverse.fasta");
		File featureFile=new File("/Volumes/BriansSSD/hela_serum_timecourse/wide_timecourse_concatenated.features.txt");
		File outputFile=new File("/Volumes/BriansSSD/hela_serum_timecourse/wide_timecourse_concatenated.peptide.txt");
		File decoyFile=new File("/Volumes/BriansSSD/hela_serum_timecourse/wide_timecourse_concatenated.peptide_decoy.txt");
		File outputProteinFile=new File("/Volumes/BriansSSD/hela_serum_timecourse/wide_timecourse_concatenated.protein.txt");
		File decoyProteinFile=new File("/Volumes/BriansSSD/hela_serum_timecourse/wide_timecourse_concatenated.protein_decoy.txt");
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
		Files.copy(is, featureFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
		
		is=getClass().getResourceAsStream("/ecoli-190209-contam_correctNL.fasta");
		File fastaFile=File.createTempFile("ecoli", ".fasta");
		fastaFile.deleteOnExit();
		Files.copy(is, fastaFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
		
		PercolatorExecutionData percolatorFiles=getPercolatorFiles(featureFile, fastaFile, SearchParameterParser.getDefaultParametersObject());


		Pair<ArrayList<PercolatorPeptide>, Float> origpair=PercolatorExecutor.executePercolatorTSV(getDefaultPercolaterVersion(), percolatorFiles, 0.01f);
		assertTrue(origpair.x.size()>0);
		assertTrue(origpair.y>0);
		
		Pair<ArrayList<PercolatorPeptide>, Float> pair=PercolatorReader.getPassingPeptidesFromTSV(percolatorFiles.getPeptideOutputFile(), 0.01f, false);
		assertEquals(origpair.x.size(), pair.x.size());
		assertEquals(origpair.y, pair.y, 0.001f);
		
		Pair<ArrayList<PercolatorPeptide>, Float> decoyPair=PercolatorReader.getPassingPeptidesFromTSV(percolatorFiles.getPeptideDecoyFile(), 0.01f, true);
		assertTrue(decoyPair.x.size()>0);
		assertTrue(decoyPair.x.size()<origpair.x.size()/99f);
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
