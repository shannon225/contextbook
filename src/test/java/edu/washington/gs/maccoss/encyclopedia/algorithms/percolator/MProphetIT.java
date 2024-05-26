package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;

import org.junit.Test;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class MProphetIT {
	public static void main(String[] args) throws Exception {
		//File featureFile=new File("/Volumes/MacOnlySSD/day8/2022_12_05_ID5_day2_tcells_16mzst_DIA_wtrap_50cm_01.mzML.features.txt");
		//File featureFile=new File("/Users/searleb/Documents/encyclopedia/small_file/bcs_2020jan16_600to603_hela_clib.dia.features.txt");
		File featureFile=new File("/Users/searleb/Documents/encyclopedia/small_file/bcs_2020jan16_hela_clib_3.dia.features.txt");
		File fastaFile=new File(featureFile.getParent(), "uniprot-9606.fasta");
		processMProphet(featureFile, fastaFile);
	}
	
	@Test
	public void testMProphet() throws Exception {
		doMProphetTest();
	}

	protected void doMProphetTest() throws IOException, InterruptedException {
		InputStream is=getClass().getResourceAsStream("/pecan.feature.txt");
		File featureFile=File.createTempFile("pecan", ".feature");
		featureFile.deleteOnExit();
		Files.copy(is, featureFile.toPath(), StandardCopyOption.REPLACE_EXISTING);

		is=getClass().getResourceAsStream("/ecoli-190209-contam_correctNL.fasta");
		File fastaFile=File.createTempFile("ecoli", ".fasta");
		fastaFile.deleteOnExit();
		Files.copy(is, fastaFile.toPath(), StandardCopyOption.REPLACE_EXISTING);

		processMProphet(featureFile, fastaFile);
	}

	private static void processMProphet(File featureFile, File fastaFile)
			throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		MProphetExecutionData percolatorFiles=getMProphetFiles(featureFile, fastaFile, SearchParameterParser.getDefaultParametersObject());

		final AminoAcidConstants aaConstants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
		final float threshold = 0.01f;
		
		MProphetResult origpair=MProphetReiter.executeMProphetTSV(percolatorFiles, threshold, aaConstants, 0);

		assertTrue(origpair.getPassingPeptides().size()>0);
		assertTrue(origpair.getPi0()>0);
		
		System.out.println("total pep: "+origpair.getPassingPeptides().size());
		System.out.println("pi_0: "+origpair.getPi0());
		
		// Check that re-reading the results gives the same data as the return from executing Percolator
		Pair<ArrayList<PercolatorPeptide>, Float> pair= PercolatorReader.getPassingPeptidesFromTSV(percolatorFiles.getPeptideOutputFile(), threshold, aaConstants, false);
		assertEquals(origpair.getPassingPeptides().size(), pair.x.size());
		assertEquals(origpair.getPi0(), pair.y, 0.001f);

		// Check for a sensible pi0
		final float pi0 = origpair.getPi0();
		assertTrue("Got invalid pi0 from MProphet (" + pi0 + ")", 0.1 < pi0 && pi0 < 0.9);

		Pair<ArrayList<PercolatorPeptide>, Float> decoyPair=PercolatorReader.getPassingPeptidesFromTSV(percolatorFiles.getPeptideDecoyFile(), threshold, aaConstants, true);
		// assert there was at least one decoy
		final int nDecoys = decoyPair.x.size();
		assertTrue(nDecoys > 0);

		// check that the decoys/targets is less than the qvalue threshold
		final int nTargets = origpair.getPassingPeptides().size();
		final float fdr = pi0 * nDecoys / (float) nTargets;
		assertTrue(
				String.format("Result didn't meet threshold! %.03f * %d / %d = %.02f >= %.02f",
						pi0,
						nDecoys,
						nTargets,
						fdr,
						threshold
				),
				fdr < threshold + 0.001f
		);
	}

	public static MProphetExecutionData getMProphetFiles(File featureFile, File fastaFile, SearchParameters parameters) throws IOException {
		File outputFile=File.createTempFile("mprophet", ".txt");
		outputFile.deleteOnExit();
		File decoyFile=File.createTempFile("mprophet", ".decoy.txt");
		decoyFile.deleteOnExit();
		MProphetExecutionData data=new MProphetExecutionData(featureFile, fastaFile, outputFile, decoyFile, parameters);
		return data;
	}
}
