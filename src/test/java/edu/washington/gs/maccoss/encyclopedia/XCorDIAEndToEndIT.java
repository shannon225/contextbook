package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.ImmutableList;
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
import java.util.List;
import java.util.Optional;

import static edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils.getResourceAsTempFile;

public class XCorDIAEndToEndIT extends AbstractEndToEndIT{
	static final String REFERENCE_SEARCH1_RESOURCE = "/edu/washington/gs/maccoss/encyclopedia/reference_data/xcorr_dia_1.elib";
	static final String REFERENCE_SEARCH2_RESOURCE = "/edu/washington/gs/maccoss/encyclopedia/reference_data/xcorr_dia_2.elib";
	static final String REFERENCE_SEARCH3_RESOURCE = "/edu/washington/gs/maccoss/encyclopedia/reference_data/xcorr_dia_3.elib";
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

		copyJobDataToResultsDirectory(jobDataA, REFERENCE_SEARCH1_RESOURCE);
		copyJobDataToResultsDirectory(jobDataB, REFERENCE_SEARCH2_RESOURCE);
		copyJobDataToResultsDirectory(jobDataC, REFERENCE_SEARCH3_RESOURCE);
	}

	@AfterClass
	public static void tearDownReports() throws Exception {
		tearDownClass();
		parameters = null;
		factory = null;
	}

	public static XCorDIAJobData makeAndDoJob(File dia) throws Exception {
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
	public List<String> getReferenceSearchResources() throws Exception {
		return ImmutableList.of(
				REFERENCE_SEARCH1_RESOURCE,
				REFERENCE_SEARCH2_RESOURCE,
				REFERENCE_SEARCH3_RESOURCE
		);
	}

	@Override
	public String getReferenceSingleQuantResource() throws Exception {
		return REFERENCE_SINGLE_QUANT_RESOURCE;
	}

	@Override
	public String getReferenceMultiResource() throws Exception {
		return REFERENCE_MULTI_RESOURCE;
	}

	@Override
	public String getReferenceMultiQuantResource() throws Exception {
		return REFERENCE_MULTI_QUANT_RESOURCE;
	}
}
