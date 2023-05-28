package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.TwoDimensionalKDE;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;

public class AlignGPFToSingleInjection {

	public static void main(String[] args) throws Exception {
		File singleInjectionForAlignment=new File("/Users/searleb/Documents/encyclopedia/alignment/2023_05_09_lcmv_pool_50kcells_16mzst_DIA_25cm_02.mzML.elib");
		File combinedGPFSaveFile=new File("/Users/searleb/Documents/encyclopedia/alignment/lcmv_pool_GPF_50kcells_combined.elib");
		
		File[] gpfs=new File[] {
			new File("/Users/searleb/Documents/encyclopedia/alignment/2023_05_09_lcmv_pool_GPF_50k_tcells_16mzst_25cm_01.mzML.elib"),
			new File("/Users/searleb/Documents/encyclopedia/alignment/2023_05_09_lcmv_pool_GPF_50k_tcells_16mzst_25cm_02.mzML.elib"),
			new File("/Users/searleb/Documents/encyclopedia/alignment/2023_05_09_lcmv_pool_GPF_50k_tcells_16mzst_25cm_03.mzML.elib"),
			new File("/Users/searleb/Documents/encyclopedia/alignment/2023_05_09_lcmv_pool_GPF_50k_tcells_16mzst_25cm_04.mzML.elib"),
			new File("/Users/searleb/Documents/encyclopedia/alignment/2023_05_09_lcmv_pool_GPF_50k_tcells_16mzst_25cm_05.mzML.elib"),
			new File("/Users/searleb/Documents/encyclopedia/alignment/2023_05_19_lcmv_pool_GPF_50k_tcells_16mzst_25cm_06.mzML.elib")
		};
		
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		
		createLibraryFromGPF(gpfs, singleInjectionForAlignment, combinedGPFSaveFile, params);
	}

	/**
	 * NOTE, THIS DOES NOT FDR CORRECT! As a result, it shouldn't be used in production
	 * @param gpfs
	 * @param singleInjectionForAlignment
	 * @param combinedGPFSaveFile
	 * @param params
	 * @throws IOException
	 * @throws SQLException
	 * @throws DataFormatException
	 */
	private static void createLibraryFromGPF(File[] gpfs, File singleInjectionForAlignment, File combinedGPFSaveFile, SearchParameters params) throws IOException, SQLException, DataFormatException {
		LibraryFile single=new LibraryFile();
		single.openFile(singleInjectionForAlignment);
		
		ArrayList<LibraryEntry> entries=single.getAllEntries(false, params.getAAConstants());
		HashMap<String, LibraryEntry> entryMap=new HashMap<>();
		for (LibraryEntry entry : entries) {
			entryMap.put(entry.getPeptideModSeq(), entry);
		}

		ArrayList<LibraryEntry> adjustedEntries=new ArrayList<>();
		for (File gpfFile : gpfs) {
			LibraryFile gpf=new LibraryFile();
			gpf.openFile(gpfFile);
			ArrayList<LibraryEntry> gpfEntries=gpf.getAllEntries(false, params.getAAConstants());
			
			ArrayList<XYPoint> rts=new ArrayList<>();
			for (LibraryEntry entry : gpfEntries) {
				LibraryEntry ref=entryMap.get(entry.getPeptideModSeq());
				if (ref!=null) {
					rts.add(new XYPoint(ref.getRetentionTime(), entry.getRetentionTime()));
				}
			}
			RetentionTimeFilter filter=RetentionTimeFilter.getFilter(rts, "Single", "GPF", TwoDimensionalKDE.HIGHER_RESOLUTION);
			filter.plot(rts, Optional.of(gpfFile));

			for (LibraryEntry entry : gpfEntries) {
				float xRT=filter.getXValue(entry.getRetentionTime());
				adjustedEntries.add(entry.updateRetentionTime(xRT));
			}
		}
		
		LibraryFile combined=new LibraryFile();
		combined.openFile();
		combined.dropIndices();
		combined.addEntries(adjustedEntries);
		combined.addProteinsFromEntries(adjustedEntries);
		combined.createIndices();
		combined.saveAsFile(combinedGPFSaveFile);
		combined.close();
	}

}
