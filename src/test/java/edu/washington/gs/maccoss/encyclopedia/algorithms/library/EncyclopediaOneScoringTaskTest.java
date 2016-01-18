package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.Encyclopedia;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.mock.MockLibrary;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;

public class EncyclopediaOneScoringTaskTest {
	public static void main(String[] args) throws Exception {

		File diaFile=new File("/Users/searleb/Documents/school/projects/pecandata/DIA_1xGFP_20x20mz_500to900_rep1.mzML");
		File libraryFile=new File("/Users/searleb/Documents/school/projects/pecandata/cptac2_human_hcd_selected.elib");
		File outputFile=new File(diaFile.getAbsolutePath()+".pecan.txt");
		File featureFile=new File(outputFile.getAbsolutePath()+".features.txt");

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		LibraryScoringFactory taskFactory=new EncyclopediaOneScoringFactory(parameters, featureFile);

		LibraryFile library=new LibraryFile();
		library.openFile(libraryFile);
		
		ArrayList<LibraryEntry> entries=library.getEntries("VDIDAPDVEVHDPDWHLK");
		MockLibrary mockLib=new MockLibrary(entries.toArray(new LibraryEntry[entries.size()]));
		
		StripeFileInterface stripefile=MzmlToDIAConverter.getFile(diaFile, parameters);
		

		Encyclopedia.runSearch(new EmptyProgressIndicator(), mockLib, stripefile, featureFile, outputFile, taskFactory);
	}

}
