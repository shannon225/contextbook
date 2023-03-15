package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.algorithms.SearchTestSupport;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import junit.framework.TestCase;

public class RepeatedFileReadingTest extends TestCase {
	private static SearchParameters PARAMETERS=SearchParameterParser.getDefaultParametersObject();
	
	public void testRead() throws Exception {
		for (int i = 0; i < 100; i++) {
			StripeFile dia=(StripeFile)StripeFileGenerator.getFile(SearchTestSupport.getNewTestSingleWindowFile(), PARAMETERS);
			LibraryFile library = (LibraryFile)BlibToLibraryConverter.getFile(SearchTestSupport.getNewTestLibraryFile());
			assertNotNull(dia);
			assertNotNull(library);
			
			library.addMetadata("key", "value");
			
			File tempFile=File.createTempFile("test", ".dlib");
			library.saveAsFile(tempFile);
			library.close();

			tempFile.delete();

			dia.addMetadata("key", "value");
			
			File tempFile2=File.createTempFile("test", ".dia");
			dia.saveAsFile(tempFile2);
			dia.close();

			tempFile.delete();
		}
	}
}
