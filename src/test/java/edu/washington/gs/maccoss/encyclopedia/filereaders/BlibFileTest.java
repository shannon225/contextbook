package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import junit.framework.TestCase;

public class BlibFileTest extends TestCase {
	public static final SearchParameters params=SearchParameterParser.getDefaultParametersObject();
	public static void main(String[] args) throws Exception {
		LibraryFile elib=new LibraryFile();
		elib.openFile(new File("/Users/searleb/Downloads/yeast/yeast_narrow_dda_library.elib"));
		ArrayList<LibraryEntry> entries=elib.getAllEntries(false, new AminoAcidConstants());
		
		SearchJobData job=new SearchJobData() {
			
			@Override
			public boolean hasBeenRun() {
				return true;
			}
			
			@Override
			public String getVersion() {
				return "Elib to Blib";
			}
			
			@Override
			public String getSearchType() {
				return "EncyclopeDIA";
			}
			
			@Override
			public PercolatorExecutionData getPercolatorFiles() {
				return null;
			}
			@Override
			public SearchJobData updateQuantFile(File f) {
				return null;
			}
			
			@Override
			public SearchParameters getParameters() {
				return params;
			}
			
			@Override
			public StripeFileInterface getDiaFileReader() {
				return null;
			}
			@Override
			public String getOriginalDiaFileName() {
				return "Dummy";
			}
			@Override
			public String getPrimaryScoreName() {
				return "primary";
			}
		};

		File blibFile=new File("/Users/searleb/Downloads/yeast/yeast_narrow_dda_library.blib");
		BlibFile blib=new BlibFile();
		blib.openFile();
		blib.setUserFile(blibFile);
		blib.dropIndices();
		blib.addLibrary(job, entries, 0, 0, 0);
		blib.createIndices();
		blib.saveFile();
		blib.close();
	}
	
	public void testCreateTables() throws IOException, SQLException {
		BlibFile file=new BlibFile();
		file.openFile();
		
		assertNotNull(file); // just getting here without errors makes sure the table structure is ok
	}
}
