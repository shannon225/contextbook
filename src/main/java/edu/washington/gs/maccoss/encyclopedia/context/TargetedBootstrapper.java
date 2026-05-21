package edu.washington.gs.maccoss.encyclopedia.context;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.WindowData;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;

public class TargetedBootstrapper {

	public static void main(String[] args) throws Throwable {
//		if (args.length!=5) {
//			Logger.errorLine("TargetedBootstrapper requires five parameters in order:");
//			Logger.logLine(" 1) Input Library (.elib file)");
//			Logger.logLine(" 2) Input raw file (.raw, .mzML, or .dia");
//			Logger.logLine(" 3) Output file path");
//			Logger.logLine(" 4) Number of peptides per assay (default = 100)");
//			Logger.logLine(" 5) Starting Seed Number (default is 3, which will use seeds 0, 1, 2, and 3, producing 4 assays total);");
		
		
//		File libraryPath = new File(args[0]);
//		File rawFile = new File(args[1]);
//		File outputPath = new File(args[2]);
//		int numberOfPeptides = Integer.parseInt(args[4]);
//		int seed = Integer.parseInt(args[5]);

		String libraryPath = "C:/Users/m334793/Documents/Library/bootstrapping_immune_cells/cd4_library.elib";
		String rawFilePath = "C:/Users/m334793/Documents/Library/bootstrapping_immune_cells/2026_01_28_EP5_CD4_PicoChipHT_30min_GPFDIA_combined00_01.dia";
		Path mapOutputPath = Paths.get("C:/Users/m334793/Documents/Library/bootstrapping_immune_cells/CD4_target_decoy_map.txt");

		Path rawFile = Paths.get(rawFilePath);
		String baseName = rawFilePath.replaceFirst("\\.dia$",  "");
		
		
		int seed = 0;
		AminoAcidConstants aaConstants = new AminoAcidConstants();
		int numberOfPeptides = 100; // number of Peptides per assay


		for (int i = 0; i <= seed; i++) {// Randomly Select Precursors, then use them to mask the .DIA file

			ArrayList<IsolationWindow> isolationWindows = selectMask(numberOfPeptides, aaConstants, i, libraryPath, mapOutputPath);
			Path outputPath = rawFile.getParent().resolve(baseName + "_masked" + i + "_assay.dia");
			Path maskedAssayOutputPath = rawFile.getParent().resolve(baseName + "_masked" + i + "_assay.txt");

			StripeFile maskedFile = writeMaskedFile(isolationWindows, i, rawFilePath, outputPath);
		    writeAssayList(isolationWindows, maskedAssayOutputPath);

			System.out.println("Complete! The masked file " + maskedFile + i + " was made.\n");
		}
	}


	// First function - Randomly Selects Precursors from a library and compiles them
	// into a list

	public static ArrayList<IsolationWindow> selectMask(int numberOfPeptides,
			AminoAcidConstants aaConstants, int i, String libraryPath, Path mapOutputPath) throws IOException, SQLException, Throwable {

		// START TIMER 1
		long startTime = System.nanoTime();
		LibraryFile library = new LibraryFile();
		File file = new File(libraryPath);
		
		ArrayList<IsolationWindow> isolationWindows = new ArrayList<>();
		
		// For mapping targets and decoys later
		HashMap<String, String> targetDecoyOriginMap = new HashMap<>();

		int randomValue = 0 + i; // Add haliburton's number to get random number
		
		HashSet<Integer> simulatedAssaySet = new HashSet<>();
		HashSet<String> sequencesSelectedForMasking = new HashSet<>();
		AminoAcidConstants constants = new AminoAcidConstants();
		SearchParameters params = PecanParameterParser.getDefaultParametersObject(); // need parameters to run smartDecoy 

		library.openFile(file);
		try {
		
			// Load all entries
			ArrayList<LibraryEntry> entries = library.getAllEntries(false, aaConstants);

			while (simulatedAssaySet.size() < numberOfPeptides) {
				randomValue = RandomGenerator.randomInt(randomValue);
				int index = Math.abs(randomValue) % entries.size();
				simulatedAssaySet.add(index);
			}
//			System.out.println("Selecting " + simulatedAssaySet.size() + " precursors for a fake assay.");

			for (Integer index : simulatedAssaySet) {
				
				// Retrieve the library entry at the random index
				LibraryEntry entry = entries.get(index);
				
				// Get the m/z, RT and sequence
				double targetMz = entry.getPrecursorMZ();
				float rtCenter = entry.getRetentionTimeInSec();
				String sequence = entry.getPeptideModSeq();
				byte charge = entry.getPrecursorCharge();


				// Calculate a RT ranges for the isolationWindows object
				float rtMin = (float) (rtCenter - (60 * 2.5));
				float rtMax = (float) (rtCenter + (60 * 2.5));
				
				// Add sequences to the isolationWindows object
				IsolationWindow window = new IsolationWindow(sequence, targetMz, charge, rtMin, rtMax, false);
				isolationWindows.add(window);
				sequencesSelectedForMasking.add(sequence);
				
				// Now let's make decoys for each target 
//				byte charge = PeptideUtils.getExpectedChargeState(sequence);
				String decoy = PeptideUtils.getSmartDecoy(sequence, charge, sequencesSelectedForMasking, params);
				String correctedDecoyMass = PeptideUtils.getCorrectedMasses(decoy, constants);
				double decoyMz = constants.getChargedMass(correctedDecoyMass, charge);
				targetDecoyOriginMap.put(sequence, decoy);   
				writeTargetDecoyMap(targetDecoyOriginMap, mapOutputPath);
				
				// Add decoys to Isolation Windows
				IsolationWindow decoyWindow = new IsolationWindow(decoy, decoyMz, charge, rtMin, rtMax, true);
				isolationWindows.add(decoyWindow);
			}
			
			library.close();
			System.out.println(isolationWindows.size() + " Precursors marked for extraction.");

		} catch (Exception e) {
			System.out.println("There was an error with selecting precursors. Check file path.");
			e.printStackTrace();
		}

		// END TIMER 1
		long endTime = System.nanoTime();
		long duration = endTime - startTime;
		System.out.println("randomlySelectPrecursors(): Time taken (ms) : " + duration / 1_000_000);

		return isolationWindows;
	}

	// Second function - Uses the IsolationWindow List to mask the raw data
	@SuppressWarnings("unused")
	public static StripeFile writeMaskedFile(ArrayList<IsolationWindow> isolationWindows, int i, String diaFilePath, Path outputPath)
			throws Throwable {

		// START TIMER 2
		long startTime = System.nanoTime();

		File rawFile = new File(diaFilePath);
		StripeFile maskedFile = new StripeFile(false);
		StripeFile rawLibraryFile = new StripeFile(false);
		File outputFile = outputPath.toFile();
		
		HashSet<Integer> addedPrecursors = new HashSet<>();
		HashSet<Integer> addedFragments = new HashSet<>();

		// System.out.println("Is the .dia file open? " + rawLibraryFile.isOpen());

		try {
			rawLibraryFile.openFile(rawFile);
			maskedFile.openFile();
			
			// Add Ranges
			HashMap<Range, WindowData> dutyCycleMap = new HashMap<>();
			System.out.println("Masking DIA file based on the selected precursors...");
			for (IsolationWindow window : isolationWindows) {	

				double windowMz = window.getTargetMz();
				float windowStartTime = window.getRtMin();
				float windowStopTime = window.getRtMax();
				boolean sqrt = false;
				double mzStart = windowMz - 1; 
				double mzStop = windowMz + 1;
				Range mzRange = new Range(mzStart, mzStop);

				ArrayList<FragmentScan> fragmentScansFromWindow = rawLibraryFile.getStripes(windowMz, windowStartTime, windowStopTime, sqrt);
				ArrayList<FragmentScan> matchingScans = new ArrayList<>();

				// Add Fragment Scans
				for (FragmentScan scan : fragmentScansFromWindow) {
					double scanMz = scan.getPrecursorMZ();
					float scanRT = scan.getScanStartTime();
					int scanIndex = scan.getSpectrumIndex();
					if (mzRange.contains(scanMz) && !addedFragments.contains(scanIndex)) {
						matchingScans.add(scan);
						addedFragments.add(scanIndex);
					} 
				}

				for (Entry<Range, WindowData> entry : rawLibraryFile.getRanges().entrySet()) { //loadRanges instead? 
					if (mzRange.contains(entry.getKey().getMiddle())) {
						dutyCycleMap.put(entry.getKey(), entry.getValue());
					}
				}
				maskedFile.setRanges(dutyCycleMap);
				maskedFile.addStripe(matchingScans);
				
				// Add Precursor Scans
				ArrayList<PrecursorScan> precursorScanFromWindow = rawLibraryFile.getPrecursors(windowStartTime, windowStopTime);
				ArrayList<PrecursorScan> matchingPrecursors = new ArrayList<>();

				for (PrecursorScan precursor : precursorScanFromWindow) {
					Range precursorRange = new Range(precursor.getIsolationWindowLower(),
							precursor.getIsolationWindowUpper());
					int spectrumIndex = precursor.getSpectrumIndex();
					if ((precursorRange.contains(mzRange) && !addedPrecursors.contains(spectrumIndex))) {
						matchingPrecursors.add(precursor);
						addedPrecursors.add(spectrumIndex);
					}
				}
				maskedFile.addPrecursor(matchingPrecursors);

			}
			maskedFile.setFileName(rawFile.getName(), null, rawFile.getAbsolutePath());
			maskedFile.addMetadata(diaFilePath, diaFilePath);
			rawLibraryFile.close();

		} catch (IOException e) {
			System.out.println("Unable to open raw file.");
			e.printStackTrace();
		}

		// END TIMER 2
		long endTime = System.nanoTime();
		long duration = endTime - startTime;

		System.out.println("maskDIAFileBasedOnIsolationWindows(): Time taken (ms) : " + duration / 1_000_000);

		maskedFile.saveAsFile(outputFile);
		System.out.println("Target mass list for the masked file  was written to " + outputPath
				+ "\n Number of added Precursor scans: " + addedPrecursors.size()
				+ "\n Number of added Fragment scans: " + addedFragments.size());
		

		return maskedFile;
	}
	
	public static void writeAssayList(ArrayList<IsolationWindow> isolationWindows, Path outputPath) throws IOException {
	    try (BufferedWriter writer = Files.newBufferedWriter(outputPath)) {
	    	writer.write("Compound\tFormula\tAdduct\tm/z\tz\tRT Time (min)\tWindow (min)\tisDecoy");
	    	writer.newLine();

	        for (IsolationWindow window : isolationWindows) {
	            String compound = window.getCompound();
	            double targetMz = window.getTargetMz();
	            byte charge = window.getCharge();
	            boolean isDecoy = window.isDecoy();

	            float rtCenterMin = ((window.getRtMin() + window.getRtMax()) / 2.0f) / 60.0f;
	            float windowMin = (window.getRtMax() - window.getRtMin()) / 60.0f;
	            
	            writer.write(compound + 
	            		"\t" +
	            		"" + 
	            		"\t" +
	            		"(no adduct)" +
	            		"\t" +
	            		targetMz +
	            		"\t" +
	            		charge +
	            		"\t" +
	            		rtCenterMin +
	            		"\t" +
	            		windowMin +
	            		"\t" +
	            		isDecoy);
	            writer.newLine();

	        } 
	    }
	}
	
	
	private static void writeTargetDecoyMap(HashMap<String, String> targetDecoyMap, Path outputPath)
	        throws IOException {

	    try (BufferedWriter writer = Files.newBufferedWriter(outputPath)) {

	        writer.write("decoySequence\ttargetSequence");
	        writer.newLine();

	        for (Entry<String, String> entry : targetDecoyMap.entrySet()) {
	            String decoySequence = entry.getKey();
	            String targetSequence = entry.getValue();

	            writer.write(decoySequence + "\t" + targetSequence);
	            writer.newLine();
	        }
	    }
	}
}
