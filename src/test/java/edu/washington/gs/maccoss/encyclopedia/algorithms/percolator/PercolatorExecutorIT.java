package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import gnu.trove.map.hash.TCharDoubleHashMap;

import org.junit.Ignore;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;

import static edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutorTest.getPercolatorFiles;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class PercolatorExecutorIT {
	@Test
	public void testPercolatorExecutorV2() throws Exception {
		doPercolatorTest(PercolatorVersion.v2p10);
	}

	@Ignore
	public void testPercolatorExecutorV3p5() throws Exception {
		doPercolatorTest(PercolatorVersion.v3p05);
	}

	@Test
	public void testPercolatorExecutorV3p1() throws Exception {
		doPercolatorTest(PercolatorVersion.v3p01);
	}

	protected void doPercolatorTest(PercolatorVersion percolatorVersion) throws IOException, InterruptedException {
		InputStream is=getClass().getResourceAsStream("/pecan.feature.txt");
		File featureFile=File.createTempFile("pecan", ".feature");
		featureFile.deleteOnExit();
		Files.copy(is, featureFile.toPath(), StandardCopyOption.REPLACE_EXISTING);

		is=getClass().getResourceAsStream("/ecoli-190209-contam_correctNL.fasta");
		File fastaFile=File.createTempFile("ecoli", ".fasta");
		fastaFile.deleteOnExit();
		Files.copy(is, fastaFile.toPath(), StandardCopyOption.REPLACE_EXISTING);

		PercolatorExecutionData percolatorFiles=getPercolatorFiles(featureFile, fastaFile, SearchParameterParser.getDefaultParametersObject());

		final AminoAcidConstants aaConstants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
		final float threshold = 0.01f;

		Pair<ArrayList<PercolatorPeptide>, Float> origpair=PercolatorExecutor.executePercolatorTSV(percolatorVersion, percolatorFiles, threshold, aaConstants);
		assertTrue(origpair.x.size()>0);
		assertTrue(origpair.y>0);

		// Check that re-reading the results gives the same data as the return from executing Percolator
		Pair<ArrayList<PercolatorPeptide>, Float> pair= PercolatorReader.getPassingPeptidesFromTSV(percolatorFiles.getPeptideOutputFile(), threshold, aaConstants, false);
		assertEquals(origpair.x.size(), pair.x.size());
		assertEquals(origpair.y, pair.y, 0.001f);

		// Check for a sensible pi0
		final float pi0 = origpair.y;
		assertTrue("Got invalid pi0 from percolator " + percolatorVersion + " (" + pi0 + ")", 0.1 < pi0 && pi0 < 0.9);

		Pair<ArrayList<PercolatorPeptide>, Float> decoyPair=PercolatorReader.getPassingPeptidesFromTSV(percolatorFiles.getPeptideDecoyFile(), threshold, aaConstants, true);
		// assert there was at least one decoy
		final int nDecoys = decoyPair.x.size();
		assertTrue(nDecoys > 0);

		// check that the decoys/targets is less than the qvalue threshold
		final int nTargets = origpair.x.size();
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
}
