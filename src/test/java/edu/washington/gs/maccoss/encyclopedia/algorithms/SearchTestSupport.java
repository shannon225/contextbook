package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;

public class SearchTestSupport {
	private static SearchParameters PARAMETERS=SearchParameterParser.getDefaultParametersObject();
	private static File SINGLE_WINDOW_FILE=null;
	private static File FULL_WINDOW_FILE=null;
	private static File LIBRARY_FILE=null;
	private static File RESULT_LIBRARY_FILE=null;
	private static File FASTA_FILE=null;
	
	public static ArrayList<FastaEntryInterface> getFastaEntries() throws Exception {
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(getTestFastaFile(), PARAMETERS);
		return entries;
	}
	
	public static LibraryFile getResultLibrary() throws Exception {
		LibraryFile library=(LibraryFile)BlibToLibraryConverter.getFile(getResultLibraryFile());
		return library;
	}
	
	public static LibraryFile getLibrary() throws Exception {
		LibraryFile library=(LibraryFile)BlibToLibraryConverter.getFile(getTestLibraryFile());
		return library;
	}
	
	public static StripeFileInterface getSingleWindowStripeFile() throws Exception {
		return StripeFileGenerator.getFile(getTestSingleWindowFile(), PARAMETERS);
	}
	
	public static StripeFileInterface getFullWindowStripeFile() throws Exception {
		return StripeFileGenerator.getFile(getTestFullWindowFile(), PARAMETERS);
	}
	
	public static SearchParameters getParameters() {
		return PARAMETERS;
	}
	
	public static File getTestSingleWindowFile() throws Exception {
		if (SINGLE_WINDOW_FILE!=null) return SINGLE_WINDOW_FILE;
		SINGLE_WINDOW_FILE = getNewTestSingleWindowFile();
		return SINGLE_WINDOW_FILE;
	}

	public static File getNewTestSingleWindowFile() throws Exception {
		return writeTempFile("/small_regression/bcs_2020jan16_600to603_hela_48to63.dia");
	}
	
	public static File getTestFullWindowFile() throws Exception {
		if (FULL_WINDOW_FILE!=null) return FULL_WINDOW_FILE;
		FULL_WINDOW_FILE = getNewTestFullWindowFile();
		return FULL_WINDOW_FILE;
	}

	public static File getNewTestFullWindowFile() throws Exception {
		return writeTempFile("/small_regression/bcs_2020jan16_hela_48p0_48p1.dia");
	}
	
	public static File getResultLibraryFile() throws Exception {
		if (RESULT_LIBRARY_FILE!=null) return RESULT_LIBRARY_FILE;
		RESULT_LIBRARY_FILE = getNewResultLibraryFile();
		return RESULT_LIBRARY_FILE;
	}
	
	public static File getTestLibraryFile() throws Exception {
		if (LIBRARY_FILE!=null) return LIBRARY_FILE;
		LIBRARY_FILE = getNewTestLibraryFile();
		return LIBRARY_FILE;
	}

	public static File getNewResultLibraryFile() throws Exception {
		return writeTempFile("/small_regression/bcs_2020jan16_600to603_hela_48to63.dia.elib");
	}

	public static File getNewTestLibraryFile() throws Exception {
		return writeTempFile("/small_regression/pan_human_library_600to603.dlib");
	}
	
	public static File getTestFastaFile() throws Exception {
		if (FASTA_FILE!=null) return FASTA_FILE;
		FASTA_FILE = getNewTestFastaFile();
		return FASTA_FILE;
	}

	public static File getNewTestFastaFile() throws Exception {
		return writeTempFile("/small_regression/pan_human_library_600to603.fasta");
	}
	
	public static File writeTempFile(String fileResourceName) throws Exception {
		File f=new File(fileResourceName);
		
		InputStream is=SearchTestSupport.class.getResourceAsStream(fileResourceName);
		return writeTempFile(is, f.getName());
	}

	/**
	 * @param is
	 * @return
	 * @throws Exception
	 */
	public static File writeTempFile(InputStream is, String filename) throws Exception {
		String tDir = System.getProperty("java.io.tmpdir");

		File output = new File(tDir, filename);
		
		output.deleteOnExit();
			Files.copy(is, output.toPath(), StandardCopyOption.REPLACE_EXISTING);
		
		return output;
	}

}
