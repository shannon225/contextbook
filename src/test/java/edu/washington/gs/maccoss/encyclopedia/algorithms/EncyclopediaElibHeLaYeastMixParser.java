package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.CoefficientOfVariationCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.LibraryReportExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.SampleCoordinate;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroupInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;

public class EncyclopediaElibHeLaYeastMixParser {
	public static String[] sampleNames=new String[] {"100.0% HeLa", "84.3% HeLa", "58.3% HeLa", "37.4% HeLa", "20.4% HeLa", "0.0% HeLa"};
	public static HashMap<String, SampleCoordinate> sampleKey=new HashMap<>();

	public static void main(String[] args) throws IOException, SQLException, DataFormatException {
		loadMap();
		
		File file=new File("/Users/searleb/Documents/chromatogram_library_manuscript/yeast_curve/hela_yeast_quants.elib");

		CoefficientOfVariationCalculator cvCalculator=new CoefficientOfVariationCalculator(sampleKey, sampleNames, 0.2f);

		LibraryFile library=new LibraryFile();
		library.openFile(file);
		
		ArrayList<ProteinGroupInterface> proteinGroups=library.getProteinGroups();
		System.out.println(proteinGroups.size()+" total protein groups");
		LibraryReportExtractor.extractMatrix(library, proteinGroups, Optional.of(cvCalculator));
	}

	public static void loadMap() {
		sampleKey.put("20170901_CalCurves_HelaCurve_N_0ug-ul_BY4741_060.mzML", new SampleCoordinate(0,0));
		sampleKey.put("20170901_CalCurves_HelaCurve_N_0ug-ul_BY4741_132.mzML", new SampleCoordinate(1,0));
		sampleKey.put("20170901_CalCurves_HelaCurve_N_0ug-ul_BY4741_197.mzML", new SampleCoordinate(2,0));
		sampleKey.put("20170901_CalCurves_HelaCurve_E_0-1ug-ul_BY4741_070.mzML", new SampleCoordinate(0,1));
		sampleKey.put("20170901_CalCurves_HelaCurve_E_0-1ug-ul_BY4741_142.mzML", new SampleCoordinate(1,1));
		sampleKey.put("20170901_CalCurves_HelaCurve_E_0-1ug-ul_BY4741_207.mzML", new SampleCoordinate(2,1));
		sampleKey.put("20170901_CalCurves_HelaCurve_D_0-3ug-ul_BY4741_072.mzML", new SampleCoordinate(0,2));
		sampleKey.put("20170901_CalCurves_HelaCurve_D_0-3ug-ul_BY4741_144.mzML", new SampleCoordinate(1,2));
		sampleKey.put("20170901_CalCurves_HelaCurve_D_0-3ug-ul_BY4741_209.mzML", new SampleCoordinate(2,2));
		sampleKey.put("20170901_CalCurves_HelaCurve_C_0-5ug-ul_BY4741_073.mzML", new SampleCoordinate(0,3));
		sampleKey.put("20170901_CalCurves_HelaCurve_C_0-5ug-ul_BY4741_145.mzML", new SampleCoordinate(1,3));
		sampleKey.put("20170901_CalCurves_HelaCurve_C_0-5ug-ul_BY4741_210.mzML", new SampleCoordinate(2,3));
		sampleKey.put("20170901_CalCurves_HelaCurve_B_0-7ug-ul_BY4741_074.mzML", new SampleCoordinate(0,4));
		sampleKey.put("20170901_CalCurves_HelaCurve_B_0-7ug-ul_BY4741_146.mzML", new SampleCoordinate(1,4));
		sampleKey.put("20170901_CalCurves_HelaCurve_B_0-7ug-ul_BY4741_211.mzML", new SampleCoordinate(2,4));
		sampleKey.put("20170901_CalCurves_HelaCurve_A_1ug-ul_BY4741_075.mzML", new SampleCoordinate(0,5));
		sampleKey.put("20170901_CalCurves_HelaCurve_A_1ug-ul_BY4741_147.mzML", new SampleCoordinate(1,5));
		sampleKey.put("20170901_CalCurves_HelaCurve_A_1ug-ul_BY4741_212.mzML", new SampleCoordinate(2,5));
	}
}
