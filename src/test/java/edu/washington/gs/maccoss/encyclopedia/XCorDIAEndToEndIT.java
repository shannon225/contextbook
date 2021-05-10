package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.ImmutableMap;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import org.junit.After;
import org.junit.Before;

import java.io.File;
import java.util.HashMap;
import java.util.Optional;

public class XCorDIAEndToEndIT extends AbstractEndToEndIT{

	PecanSearchParameters parameters;
	XCorDIAOneScoringFactory factory;

	@Before
	public void setUp() throws Exception {
		parameters = PecanParameterParser.parseParameters(new HashMap<>(ImmutableMap.of(
				"-percolatorThreshold", "0.25",
				"-percolatorTrainingFDR","0.25",
				"-percolatorProteinThreshold","0.25")));
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
