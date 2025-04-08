package edu.washington.gs.maccoss.encyclopedia.filereaders.spectrumprocessors;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Map;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.DIAProcessor;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.WindowData;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import junit.framework.TestCase;

public class WindowDownsamplerTest extends TestCase {
	public static void main(String[] args) throws Exception {
		File startFile=new File("/Users/searleb/Documents/encyclopedia/small_file/bcs_2020jan16_hela_clib_3.dia");
		File resultFile=new File("/Users/searleb/Documents/encyclopedia/small_file/asIf4Windows.dia");
		
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		
		ArrayList<Range> downsampledRanges=new ArrayList<>();
		downsampledRanges.add(new Range(598.522827148437f,	626.53466796875f));
		downsampledRanges.add(new Range(626.53466796875f,	650.545593261719f));
		downsampledRanges.add(new Range(650.545593261719f,	676.557434082031f));
		downsampledRanges.add(new Range(676.557434082031f,	702.568420410156f));
		
		WindowDownsampler downsampler=new WindowDownsampler(downsampledRanges, parameters.getFragmentTolerance());
		DIAProcessor processor=new DIAProcessor(downsampler, parameters);

		StripeFileInterface startStripeFile = StripeFileGenerator.getFile(startFile, parameters);
		
		processor.processStripeFile(new EmptyProgressIndicator(), startStripeFile, resultFile, false);
	}
	
	public void testDownsample() throws Exception {
		final Path startFile = getResourceAsTempFile("WindowDownsamplerTest", ".dia", "/small_regression/bcs_2020jan16_hela_48p0_48p1.dia");
		final File resultFile=File.createTempFile("WindowDownsamplerTestResult", ".dia");
		resultFile.deleteOnExit();
		
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		
		ArrayList<Range> downsampledRanges=new ArrayList<>();
		downsampledRanges.add(new Range(598.522827148437f,	622.532897949219f));
		downsampledRanges.add(new Range(622.532836914062f,	646.543823242188f));
		downsampledRanges.add(new Range(646.543762207031f,	670.5546875f));
		downsampledRanges.add(new Range(670.5546875f,	694.565612792969f));
		
		WindowDownsampler downsampler=new WindowDownsampler(downsampledRanges, parameters.getFragmentTolerance());
		DIAProcessor processor=new DIAProcessor(downsampler, parameters);

		StripeFileInterface startStripeFile = StripeFileGenerator.getFile(startFile.toFile(), parameters);
		
		processor.processStripeFile(new EmptyProgressIndicator(), startStripeFile, resultFile, false);
		
		StripeFileInterface file = StripeFileGenerator.getFile(resultFile, parameters);
		Map<Range, WindowData> ranges=file.getRanges();
		assertEquals(ranges.size(), downsampledRanges.size());
		for (Entry<Range, WindowData> entry : ranges.entrySet()) {
			assertTrue(downsampledRanges.indexOf(entry.getKey())>=0);
		}
		
		assertEquals(3, file.getPrecursors(0, Float.MAX_VALUE).size());
		ArrayList<FragmentScan> stripes = file.getStripes(new Range(0, 99999), 0.0f, Float.MAX_VALUE, false);
	
		for (FragmentScan fragmentScan : stripes) {
			System.out.println(fragmentScan.getScanStartTime()+"\t"+fragmentScan.getRange());
		}
		assertEquals(12, stripes.size());
	}
	
	protected static Path getResourceAsTempFile(String name, String suffix, String resource) throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(new WindowDownsamplerTest().getClass(), resource, null, name, suffix);
	}
}
