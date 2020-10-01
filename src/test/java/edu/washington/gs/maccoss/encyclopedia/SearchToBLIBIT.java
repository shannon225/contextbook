package edu.washington.gs.maccoss.encyclopedia;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.junit.Assert.*;

public class SearchToBLIBIT {
	private final SearchParameters searchParameters = SearchParameterParser.getDefaultParametersObject();
	private final ProgressIndicator progress = new EmptyProgressIndicator();
	private Path tempDir;

	@Before
	public void setUp() throws Exception {
		tempDir = Files.createTempDirectory("SearchToBLIBIT_");
		FileUtils.forceDeleteOnExit(tempDir.toFile());
	}

	@After
	public void tearDown() throws Exception {
		if (null != tempDir) {
			FileUtils.deleteDirectory(tempDir.toFile());
		}
	}

	@Test
	public void testConvertMultiSampleLibrary() throws Exception {
		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".elib");
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final List<SearchJobData> jobData = Stream.of("") //TODO: test data
				.map(name -> new EncyclopediaJobData(
						null,
						new PercolatorExecutionData(
							//TODO
						),
						searchParameters,
						"SearchToBLIBIT",
						null,
						new EncyclopediaOneScoringFactory(searchParameters)
				))
				.collect(Collectors.toList());

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				false,
				true
		);

		//TODO: assert counts of passing peptides
		//TODO: assert percolator version, pi0
	}

	//TODO: test multi-sample combination (blib)
	//TODO: test multi-sample combination (quant)
	//TODO: test single-sample combination (elib)
	//TODO: test single-sample combination (blib)
	//TODO: test single-sample combination (quant)
}