package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrerInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import org.junit.Assume;
import org.junit.Test;

public class SearchToBLIBTest {
	@Test
	public void testAlignmentRoundTrip() throws Exception {
		final SearchParameters parameters = SearchParameterParser.getDefaultParametersObject();
		final SearchJobData job = mockJob(parameters);
		final PeakLocationInferrerInterface inferrer = mockInferrer(job);

		final LibraryFile file = new LibraryFile();
		try {
			file.openFile();

			file.addRtAlignment(job, inferrer);

			final PeakLocationInferrerInterface read = SearchToBLIB.readInferrer(file, ImmutableList.of(job), parameters);

			assertSameRtWarping(inferrer, read);
		} finally {
			file.close();
		}
	}

	private SearchJobData mockJob(SearchParameters parameters) {
		return null; //TODO
	}

	private PeakLocationInferrerInterface mockInferrer(SearchJobData... jobs) {
		return null; //TODO
	}

	private void assertSameRtWarping(PeakLocationInferrerInterface expected, PeakLocationInferrerInterface actual) {
		Assume.assumeTrue("TODO: implement assertions", false); //TODO
	}
}
