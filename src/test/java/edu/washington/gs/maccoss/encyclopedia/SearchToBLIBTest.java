package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrerInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface.AlignmentDataPoint;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.WindowData;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.QuantitativeDIAData;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.*;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.zip.DataFormatException;

import static org.junit.Assert.assertArrayEquals;

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

			assertSameRtWarping(inferrer, read, job);
		} finally {
			file.close();
		}
	}

	private static void assertSameRtWarping(PeakLocationInferrerInterface expected, PeakLocationInferrerInterface actual, SearchJobData job) {
		final List<AlignmentDataPoint> sortExpected = Lists.newArrayList(expected.getAlignmentData(job));
		final List<AlignmentDataPoint> sortActual = Lists.newArrayList(actual.getAlignmentData(job));

		Collections.sort(sortExpected, Comparator.comparing(AlignmentDataPoint::getLibrary));
		Collections.sort(sortActual, Comparator.comparing(AlignmentDataPoint::getLibrary));

		for (ToDoubleFunction<AlignmentDataPoint> fn : ImmutableList.<ToDoubleFunction<AlignmentDataPoint>>of(
				AlignmentDataPoint::getLibrary,
				AlignmentDataPoint::getActual,
				AlignmentDataPoint::getPredictedActual,
				AlignmentDataPoint::getProbability
		)) {
			assertArrayEquals(
					sortExpected.stream().mapToDouble(fn).toArray(),
					sortActual.stream().mapToDouble(fn).toArray(),
					0.0001
			);
		}
	}

	private static SearchJobData mockJob(SearchParameters parameters) {
		return new QuantitativeSearchJobData(null, null, null, parameters, "test") {
			@Override
			public File getResultLibrary() {
				return null;
			}

			@Override
			public String getSearchType() {
				return "test";
			}

			@Override
			public StripeFileInterface getDiaFileReader() {
				return new StripeFileInterface() {
					@Override
					public Map<Range, WindowData> getRanges() {
						throw new UnsupportedOperationException("Not mocked");
					}

					@Override
					public void openFile(File userFile) throws IOException, SQLException {

					}

					@Override
					public ArrayList<PrecursorScan> getPrecursors(float minRT, float maxRT) throws IOException, SQLException, DataFormatException {
						throw new UnsupportedOperationException("Not mocked");
					}

					@Override
					public ArrayList<FragmentScan> getStripes(double targetMz, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException {
						throw new UnsupportedOperationException("Not mocked");
					}

					@Override
					public ArrayList<FragmentScan> getStripes(Range targetMzRange, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException {
						throw new UnsupportedOperationException("Not mocked");
					}

					@Override
					public float getTIC() throws IOException, SQLException {
						throw new UnsupportedOperationException("Not mocked");
					}

					@Override
					public float getGradientLength() throws IOException, SQLException {
						throw new UnsupportedOperationException("Not mocked");
					}

					@Override
					public void close() {

					}

					@Override
					public boolean isOpen() {
						throw new UnsupportedOperationException("Not mocked");
					}

					@Override
					public File getFile() {
						throw new UnsupportedOperationException("Not mocked");
					}

					@Override
					public String getOriginalFileName() {
						return "test_file";
					}
				};
			}
		};
	}

	private static PeakLocationInferrerInterface mockInferrer(SearchJobData... jobs) {
		return new PeakLocationInferrerInterface() {
			private final Map<SearchJobData, List<AlignmentDataPoint>> dataMap;
			{
				dataMap = Maps.newHashMap();
				for (SearchJobData job : jobs) {
					dataMap.put(job, mockAlignmentData());
				}
			}

			@Override
			public Optional<QuantitativeDIAData> getQuantitativeData(TransitionRefinementData data) {
				throw new UnsupportedOperationException("Not mocked");
			}

			@Override
			public double[] getTopNBestIons(String peptideModSeq, byte precursorCharge) {
				throw new UnsupportedOperationException("Not mocked");
			}

			@Override
			public float getPreciseRTInSec(SearchJobData job, String peptideModSeq, float detectedRTInSec) {
				throw new UnsupportedOperationException("Not mocked");
			}

			@Override
			public float getWarpedRTInSec(SearchJobData job, String peptideModSeq) {
				throw new UnsupportedOperationException("Not mocked");
			}

			@Override
			public List<AlignmentDataPoint> getAlignmentData(SearchJobData job) {
				return dataMap.get(job);
			}
		};
	}

	private static List<AlignmentDataPoint> mockAlignmentData() {
		return Stream.generate(SearchToBLIBTest::mockAlignmentPoint)
				.limit(250)
				.collect(Collectors.toList());
	}

	private static final Random random = new Random();
	private static AlignmentDataPoint mockAlignmentPoint() {
		final float lib = 100 * random.nextFloat() + 1;
		final float pred = 5 * lib - 100;
		final float delta = (1 - random.nextFloat()) * 5;
		final float actual = pred + delta;
		final float prob = 1f;

		return AlignmentDataPoint.of(
				lib,
				actual,
				pred,
				delta,
				prob,
				false,
				"fakePep"
		);
	}
}
