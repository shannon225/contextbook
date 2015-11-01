package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.DotProduct;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Swath;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SwathFile;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class SwathExtractor {
	public static void main(String[] args) throws IOException, SQLException, DataFormatException {
		File f=new File("/Users/searleb/Documents/school/projects/mzml/q06051_rl_MCF7_IMAC_GpX_6.dia");
		SwathFile swathfile=new SwathFile();
		swathfile.openFile(f);
		
		File lf=new File("/Users/searleb/Documents/school/projects/qe_phospho.elib");
		LibraryFile libraryFile=new LibraryFile();
		libraryFile.openFile(lf);
		
		MassTolerance tolerance = new MassTolerance(10);
		DotProduct scorer=new DotProduct(tolerance);
		for (Range range : swathfile.getRanges()) {
			for (LibraryEntry entry : libraryFile.getEntries(range.getStart(), range.getStop())) {
				ArrayList<Swath> swaths=swathfile.getSwaths(entry.getPrecursorMZ(), -Float.MAX_VALUE, Float.MAX_VALUE);
				
				EValueCalculator calculator = scorePeptide(scorer, entry, swaths);
				EValueCalculator reverseCalculator = scorePeptide(scorer, entry.getReverse(tolerance), swaths);
				System.out.println(calculator.getNegLog10EValue()+"\t"+reverseCalculator.getNegLog10EValue()+"\t"+entry.getPeptideModSeq());
			}
		}
	}

	private static EValueCalculator scorePeptide(DotProduct scorer,
			LibraryEntry entry, ArrayList<Swath> swaths) {
		TFloatFloatHashMap scoreMap=new TFloatFloatHashMap();
		for (Swath swath : swaths) {
			float score=scorer.score(entry, swath);
			scoreMap.put(swath.getScanStartTime(), score);
		}
		EValueCalculator calculator=new EValueCalculator(scoreMap);
		return calculator;
	}

}
