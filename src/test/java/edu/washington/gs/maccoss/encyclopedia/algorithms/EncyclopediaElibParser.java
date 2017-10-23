package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.CoefficientOfVariationCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.LibraryReportExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.SampleCoordinate;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;

public class EncyclopediaElibParser {
	public static String[] sampleNames=new String[] {"0 hr", "2 hr", "4 hr", "8 hr", "16 hr", "24 hr"};
	public static HashMap<String, SampleCoordinate> sampleKey=new HashMap<>();

	public static void main(String[] args) throws IOException, SQLException, DataFormatException {
		loadMap();
		
		File file=new File("/Volumes/BriansSSD/hela_serum_timecourse/combined.elib");

		CoefficientOfVariationCalculator cvCalculator=new CoefficientOfVariationCalculator(sampleKey, sampleNames, 0.2f);

		LibraryFile library=new LibraryFile();
		library.openFile(file);
		
		LibraryReportExtractor.extractMatrix(library, library.getProteinGroups(), Optional.of(cvCalculator));
	}

	public static void loadMap() {
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_1a.mzML", new SampleCoordinate(0, 0));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_1b.mzML", new SampleCoordinate(0, 1));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_1c.mzML", new SampleCoordinate(0, 2));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_1d.mzML", new SampleCoordinate(0, 3));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_1e.mzML", new SampleCoordinate(0, 4));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_1f.mzML", new SampleCoordinate(0, 5));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_2a.mzML", new SampleCoordinate(1, 0));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_2b.mzML", new SampleCoordinate(1, 1));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_2c.mzML", new SampleCoordinate(1, 2));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_2d.mzML", new SampleCoordinate(1, 3));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_2e.mzML", new SampleCoordinate(1, 4));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_2f.mzML", new SampleCoordinate(1, 5));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_3a.mzML", new SampleCoordinate(2, 0));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_3b.mzML", new SampleCoordinate(2, 1));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_3c.mzML", new SampleCoordinate(2, 2));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_3d.mzML", new SampleCoordinate(2, 3));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_3e.mzML", new SampleCoordinate(2, 4));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_3f.mzML", new SampleCoordinate(2, 5));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_4a.mzML", new SampleCoordinate(3, 0));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_4b.mzML", new SampleCoordinate(3, 1));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_4c.mzML", new SampleCoordinate(3, 2));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_4d.mzML", new SampleCoordinate(3, 3));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_4e.mzML", new SampleCoordinate(3, 4));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_4f.mzML", new SampleCoordinate(3, 5));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_5a.mzML", new SampleCoordinate(4, 0));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_5b.mzML", new SampleCoordinate(4, 1));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_5c.mzML", new SampleCoordinate(4, 2));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_5d.mzML", new SampleCoordinate(4, 3));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_5e.mzML", new SampleCoordinate(4, 4));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_5f.mzML", new SampleCoordinate(4, 5));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_6a.mzML", new SampleCoordinate(5, 0));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_6b.mzML", new SampleCoordinate(5, 1));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_6c.mzML", new SampleCoordinate(5, 2));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_6d.mzML", new SampleCoordinate(5, 3));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_6e.mzML", new SampleCoordinate(5, 4));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_6f.mzML", new SampleCoordinate(5, 5));
	}
}
