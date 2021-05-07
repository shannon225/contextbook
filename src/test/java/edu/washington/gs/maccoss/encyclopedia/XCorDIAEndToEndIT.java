package edu.washington.gs.maccoss.encyclopedia;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import org.junit.After;
import org.junit.Before;

import java.io.File;
import java.util.Optional;

import static edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils.getResourceAsTempFile;

public class XCorDIAEndToEndIT extends AbstractEndToEndIT{

	PecanSearchParameters parameters;
	XCorDIAOneScoringFactory factory;

	@Before
	public void setUp() throws Exception {
		PecanSearchParameters defaultParameters = PecanParameterParser.getDefaultParametersObject();
		//we need to write out the whole constructor just to change one parameter :w
		parameters = new PecanSearchParameters(defaultParameters.getAAConstants(),
				defaultParameters.getFragType(),
				defaultParameters.getPrecursorTolerance(),
				defaultParameters.getPrecursorOffsetPPM(),
				defaultParameters.getPrecursorIsolationMargin(),
				defaultParameters.getFragmentTolerance(),
				defaultParameters.getFragmentOffsetPPM(),
				defaultParameters.getEnzyme(),
				defaultParameters.getExpectedPeakWidth(),
				defaultParameters.getMinPeptideLength(),
				defaultParameters.getMaxPeptideLength(),
				defaultParameters.getMaxMissedCleavages(),
				defaultParameters.getMinCharge(),
				defaultParameters.getMaxCharge(),
				defaultParameters.getNumberOfReportedPeaks(),
				defaultParameters.isAddDecoysToBackgound(),
				defaultParameters.isDontRunDecoys(),
				0.25f,
				0.25f,
				defaultParameters.getPercolatorVersionNumber(),
				defaultParameters.getPercolatorTrainingSetSize(),
				0.25f,
				defaultParameters.getAlpha(),
				defaultParameters.getBeta(),
				defaultParameters.getDataAcquisitionType(),
				defaultParameters.getNumberOfThreadsUsed(),
				defaultParameters.getTargetWindowCenter(),
				defaultParameters.getPrecursorWindowSize(),
				defaultParameters.getNumberOfQuantitativePeaks(),
				defaultParameters.getMinNumOfQuantitativePeaks(),
				defaultParameters.getTopNTargetsUsed(),
				defaultParameters.getMinIntensity(),
				defaultParameters.isQuantifySameFragmentsAcrossSamples(),
				defaultParameters.isVerifyModificationIons(),
				defaultParameters.isRequireVariableMods(),
				defaultParameters.isFilterPeaklists(),
				defaultParameters.isDoNotUseGlobalFDR(),
				defaultParameters.isEnableAdvancedOptions()
		);
		factory = new XCorDIAOneScoringFactory(parameters);
		super.setUp();
	}

	@After
	public void tearDown() throws Exception {
		super.tearDown();
		parameters = null;
		factory = null;
	}

	@Override
	public SearchJobData makeAndDoJob(File dia) throws Exception {
		XCorDIAJobData jobData=new XCorDIAJobData(Optional.empty(), Optional.empty(), dia, fastaFile, new File(dia.getAbsolutePath()+XCorDIAJobData.OUTPUT_FILE_SUFFIX), factory);
		XCorDIA.runPie(new EmptyProgressIndicator(), jobData);
		return jobData;
	}

	@Override
	public int getPeptideFloor() {
		return 50;
	}

	@Override
	public int getProteinFloor() {
		return 50;
	}

	@Override
	public LibraryFile[] getReferenceSearches() throws Exception {
		File referenceSearch1File = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/current_reports/xcorr_dia_1.elib", tempDir, "reference1_", ".elib").toFile();
		LibraryFile referenceSearch1 = new LibraryFile() {{openFile(referenceSearch1File);}};

		File referenceSearch2File = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/current_reports/xcorr_dia_2.elib", tempDir, "reference2_", ".elib").toFile();
		LibraryFile referenceSearch2 = new LibraryFile() {{openFile(referenceSearch2File);}};

		File referenceSearch3File = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/current_reports/xcorr_dia_3.elib", tempDir, "reference3_", ".elib").toFile();
		LibraryFile referenceSearch3 = new LibraryFile() {{openFile(referenceSearch3File);}};

		return new LibraryFile[] {referenceSearch1,referenceSearch2,referenceSearch3};
	}

	@Override
	public LibraryFile getReferenceSingleQuant() throws Exception {
		File referenceSingleQuant = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/current_reports/xcorr_single.elib", tempDir, "single_", ".elib").toFile();
		return new LibraryFile() {{openFile(referenceSingleQuant);}};
	}

	@Override
	public LibraryFile getReferenceMulti() throws Exception {
		File referenceMulti = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/current_reports/xcorr_multi.elib", tempDir, "multi_", ".elib").toFile();
		return new LibraryFile() {{openFile(referenceMulti);}};
	}

	@Override
	public LibraryFile getReferenceMultiQuant() throws Exception {
		File referenceMultiQuant = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/current_reports/xcorr_multi_quant.elib", tempDir, "multi_quant_", ".elib").toFile();
		return new LibraryFile() {{openFile(referenceMultiQuant);}};
	}

}
