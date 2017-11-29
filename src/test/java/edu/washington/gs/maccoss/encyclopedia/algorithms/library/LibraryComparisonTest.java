package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Correlation;

public class LibraryComparisonTest {
	public static void main(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		
		//File libraryFile=new File("/Users/searleb/Documents/school/encyclopedia_manuscript/HeLa.elib");
		File libraryFile=new File("/Users/searleb/Documents/school/encyclopedia_manuscript/22oct2017_hela_serum_timecourse_narrow_library.elib");
		LibraryFile file=new LibraryFile();
		file.openFile(libraryFile);
		ArrayList<LibraryEntry> entries=file.getEntries(new Range(0, 10000), false);
		HashMap<String, LibraryEntry> ddaMap=new HashMap<>();
		for (LibraryEntry entry : entries) {
			ddaMap.put(getKey(entry), AnnotatedLibraryEntry.getAnnotationsOnly(entry, parameters));
		}
		file.close();
		
		//libraryFile=new File("/Users/searleb/Documents/school/encyclopedia_manuscript/dda_lib_23aug2017_hela_serum_timecourse_pool_wide_001_170829031834.mzML.elib");
		libraryFile=new File("/Users/searleb/Documents/school/encyclopedia_manuscript/correct_lib_23aug2017_hela_serum_timecourse_pool_wide_001_170829031834.mzML.elib");
		file=new LibraryFile();
		file.openFile(libraryFile);
		entries=file.getEntries(new Range(0, 10000), false);
		
		MassTolerance tolerance=new MassTolerance(10.0);
		System.out.println("Pearson Correlation Coefficient\tCharge\tPearson Correlation Coefficient");
		for (LibraryEntry dia : entries) {
			LibraryEntry dda=ddaMap.get(getKey(dia));
			if (dda==null) System.out.println("MISSING "+getKey(dia));
			
			float correlation=(float)Correlation.getPearsons(dia, dda, tolerance);
			System.out.println(correlation+"\t"+dia.getPrecursorCharge()+"\t"+(Math.round(correlation*50f)/50f));
		}
	}
	
	private static String getKey(LibraryEntry e) {
		// avoids modification issues
		return e.getPeptideSeq()+"_"+Math.round(e.getPrecursorMZ());
	}
}
