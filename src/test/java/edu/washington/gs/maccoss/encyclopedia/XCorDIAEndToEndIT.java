package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.ImmutableMap;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import org.junit.AfterClass;
import org.junit.BeforeClass;

import java.io.File;
import java.util.HashMap;
import java.util.Optional;

import static edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils.getResourceAsTempFile;

public class XCorDIAEndToEndIT extends AbstractEndToEndIT{
	static final String REFERENCE_SEARCH1_RESOURCE = "/edu/washington/gs/maccoss/encyclopedia/reference_data/xcorr_dia_1.elib";
	static final String REFERENCE_SEARCH2_RESOURCE = "/edu/washington/gs/maccoss/encyclopedia/reference_data/xcorr_dia_2.elib";
	static final String REFERENCE_SEARCH3_ELIB = "/edu/washington/gs/maccoss/encyclopedia/reference_data/xcorr_dia_3.elib";
	static final String REFERENCE_SINGLE_QUANT_RESOURCE = "/edu/washington/gs/maccoss/encyclopedia/reference_data/xcorr_single.elib";
	static final String REFERENCE_MULTI_RESOURCE = "/edu/washington/gs/maccoss/encyclopedia/reference_data/xcorr_multi.elib";
	static final String REFERENCE_MULTI_QUANT_RESOURCE = "/edu/washington/gs/maccoss/encyclopedia/reference_data/xcorr_multi_quant.elib";

	static PecanSearchParameters parameters;
	static XCorDIAOneScoringFactory factory;

	/**
	 * These values are low estimates for the number of peptides
	 * and protein groups identified by Xcordia with 0.25 training
	 * set threshold.
	 */
	static int PEPTIDE_FLOOR = 400;
	static int PROTEIN_FLOOR = 400;

	@BeforeClass
	public static void buildReports() throws Exception {
		parameters = PecanParameterParser.parseParameters(new HashMap<>(ImmutableMap.of(
				"-percolatorThreshold", "0.25",
				"-percolatorTrainingFDR","0.25",
				"-percolatorProteinThreshold","0.25")));
		factory = new XCorDIAOneScoringFactory(parameters);
		setUpClass();
		jobDataA = makeAndDoJob(diaFile);
		jobDataB = makeAndDoJob(diaFile2);
		jobDataC = makeAndDoJob(diaFile3);
	}

	@AfterClass
	public static void tearDownReports() throws Exception {
		tearDownClass();
		parameters = null;
		factory = null;
	}

	public static SearchJobData makeAndDoJob(File dia) throws Exception {
		XCorDIAJobData jobData=new XCorDIAJobData(Optional.empty(), Optional.empty(), dia, fastaFile, new File(dia.getAbsolutePath()+XCorDIAJobData.OUTPUT_FILE_SUFFIX), factory);
		XCorDIA.runPie(new EmptyProgressIndicator(), jobData);
		return jobData;
	}

	@Override
	public int getPeptideFloor() {
		return PEPTIDE_FLOOR;
	}

	@Override
	public int getProteinFloor() {
		return PROTEIN_FLOOR;
	}

	@Override
	public LibraryFile[] getReferenceSearches() throws Exception {
		final File referenceSearch1File = getResourceAsTempFile(getClass(), REFERENCE_SEARCH1_RESOURCE, tempDir, "reference1_", ".elib").toFile();
		final LibraryFile referenceSearch1 = new LibraryFile();
		referenceSearch1.openFile(referenceSearch1File);

		final File referenceSearch2File = getResourceAsTempFile(getClass(), REFERENCE_SEARCH2_RESOURCE, tempDir, "reference2_", ".elib").toFile();
		final LibraryFile referenceSearch2 = new LibraryFile();
		referenceSearch2.openFile(referenceSearch2File);

		final File referenceSearch3File = getResourceAsTempFile(getClass(), REFERENCE_SEARCH3_ELIB, tempDir, "reference3_", ".elib").toFile();
		final LibraryFile referenceSearch3 = new LibraryFile();
		referenceSearch3.openFile(referenceSearch3File);

		return new LibraryFile[] {referenceSearch1,referenceSearch2,referenceSearch3};
	}

	@Override
	public LibraryFile getReferenceSingleQuant() throws Exception {
		final File referenceSingleQuant = getResourceAsTempFile(getClass(), REFERENCE_SINGLE_QUANT_RESOURCE, tempDir, "single_", ".elib").toFile();
		final LibraryFile libraryFile = new LibraryFile();
		libraryFile.openFile(referenceSingleQuant);
		return libraryFile;
	}

	@Override
	public LibraryFile getReferenceMulti() throws Exception {
		final File referenceMulti = getResourceAsTempFile(getClass(), REFERENCE_MULTI_RESOURCE, tempDir, "multi_", ".elib").toFile();
		final LibraryFile libraryFile = new LibraryFile();
		libraryFile.openFile(referenceMulti);
		return libraryFile;
	}

	@Override
	public LibraryFile getReferenceMultiQuant() throws Exception {
		final File referenceMultiQuant = getResourceAsTempFile(getClass(), REFERENCE_MULTI_QUANT_RESOURCE, tempDir, "multi_quant_", ".elib").toFile();
		final LibraryFile libraryFile = new LibraryFile();
		libraryFile.openFile(referenceMultiQuant);
		return libraryFile;
	}

}
