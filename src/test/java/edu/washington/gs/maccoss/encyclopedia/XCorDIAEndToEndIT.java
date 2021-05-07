package edu.washington.gs.maccoss.encyclopedia;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import org.junit.After;
import org.junit.Before;

import java.io.File;
import java.util.Optional;

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
				0.85f,
				0.85f,
				defaultParameters.getPercolatorVersionNumber(),
				defaultParameters.getPercolatorTrainingSetSize(),
				defaultParameters.getPercolatorTrainingSetThreshold(),
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

}
