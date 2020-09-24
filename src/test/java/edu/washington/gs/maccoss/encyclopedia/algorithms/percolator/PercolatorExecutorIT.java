package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import gnu.trove.map.hash.TCharDoubleHashMap;
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
		doPercolatorTest(2);
	}

	@Test
	public void testPercolatorExecutorV3() throws Exception {
		doPercolatorTest(3);
	}

	protected void doPercolatorTest(int percolatorVersion) throws IOException, InterruptedException {
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

		Pair<ArrayList<PercolatorPeptide>, Float> origpair=PercolatorExecutor.executePercolatorTSV(percolatorVersion, percolatorFiles, 0.01f, aaConstants);
		assertTrue(origpair.x.size()>0);
		assertTrue(origpair.y>0);

		Pair<ArrayList<PercolatorPeptide>, Float> pair= PercolatorReader.getPassingPeptidesFromTSV(percolatorFiles.getPeptideOutputFile(), 0.01f, aaConstants, false);
		assertEquals(origpair.x.size(), pair.x.size());
		assertEquals(origpair.y, pair.y, 0.001f);

		Pair<ArrayList<PercolatorPeptide>, Float> decoyPair=PercolatorReader.getPassingPeptidesFromTSV(percolatorFiles.getPeptideDecoyFile(), 0.01f, aaConstants, true);
		assertTrue(decoyPair.x.size()>0);
		assertTrue(decoyPair.x.size()<origpair.x.size()/99f);
	}
}
