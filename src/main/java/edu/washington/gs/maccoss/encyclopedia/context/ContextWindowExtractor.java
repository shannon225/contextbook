package edu.washington.gs.maccoss.encyclopedia.context;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.WindowData;

public class ContextWindowExtractor {

	public static void main(String[] args) throws IOException, SQLException {
		String targetAssayPath = "C:/Users/m334793/Documents/Library/assay7.csv";
		Path rawFilePath = Paths.get("C:/Users/m334793/Documents/Library/masked1_cd14_combined.dia");
		
		String rawFileName = rawFilePath.getFileName().toString();
		String baseName = rawFileName.replaceFirst("\\.dia$",  "");
		
		Path outputPath = rawFilePath.getParent().resolve(baseName + "_context.dia");

		StripeFile maskedFile = maskByTargetWindows(rawFilePath, outputPath, targetAssayPath);
		System.out.println("Masked File was created at " + outputPath + "\n" + maskedFile);
	}

	public static StripeFile maskByTargetWindows(Path rawFilePath, Path outputPath, String scheduledAssayPath)
			throws IOException, SQLException {
		long startTime = System.nanoTime();

		// Make a file output from the argument 'output path'
		File outputFile = outputPath.toFile();

		// And other files
		File rawFile = rawFilePath.toFile();
		StripeFile targetFile = new StripeFile();
		StripeFile maskedFile = new StripeFile(false);

		try {

			// Define ArrayList<> of IsolationWindows 
			// If this works, you should see a list of windows being added to scheduledWindows
			ArrayList<IsolationWindow> scheduledWindows = IsolationWindowReader.parseMassList(scheduledAssayPath);

			targetFile.openFile(rawFile);
			maskedFile.openFile();

			// Add ranges to the masked file 
			HashMap<Range, WindowData> dutyCycleMap = new HashMap<>();

			HashSet<Integer> addedFragments = new HashSet<>();
			HashSet<Double> addedPrecursors = new HashSet<>();

			for (IsolationWindow window : scheduledWindows) {
				double extractionWidthMz = 0.7;
				double windowMz = window.getTargetMz();
				double mzStart = windowMz - (extractionWidthMz / 2);
				double mzStop = windowMz + (extractionWidthMz / 2);
				Range mzRange = new Range(mzStart, mzStop);
				float windowRtMin = window.getRtMin();
				float windowRtMax = window.getRtMax();

				ArrayList<FragmentScan> fragmentScansInWindow = targetFile.getStripes(mzRange, windowRtMin,
						windowRtMax, false);
				ArrayList<FragmentScan> matchingScans = new ArrayList<>();

				// Find fragment scans that match the isolation window range for m/z 
				for (FragmentScan scan : fragmentScansInWindow) {
					double scanMz = scan.getPrecursorMZ();
					addedPrecursors.add(scanMz);
					int scanIndex = scan.getSpectrumIndex();
					if (mzRange.contains(scanMz) && !addedFragments.contains(scanIndex)) {
						matchingScans.add(scan);
						addedFragments.add(scanIndex);
						// System.out.println("Scan index added " + scanIndex);
					}
				}
				maskedFile.addStripe(matchingScans);

				// Find precursors
//				ArrayList<PrecursorScan> precursorScanFromWindow = targetFile.getPrecursors(windowRtMin, windowRtMax);
				ArrayList<PrecursorScan> matchingPrecursors = new ArrayList<>();
				ArrayList<PrecursorScan> precursorScansFromWindow = targetFile.getPrecursors(windowRtMin, windowRtMax);

				for (PrecursorScan precursor : precursorScansFromWindow) {
					Range precursorRange = new Range(precursor.getScanStartTime(),
							precursor.getIsolationWindowUpper());
					double scanIndex = precursor.getSpectrumIndex();
					double[] fraction = precursor.getMassArray();
					System.out.println("Fraction: " + fraction);
					if ((precursorRange.contains(mzRange) && !addedPrecursors.contains(scanIndex))) {
						matchingPrecursors.add(precursor);
						addedPrecursors.add(scanIndex);
					}
				}
				maskedFile.addPrecursor(matchingPrecursors);

				// Load Ranges from the raw file.
				for (Entry<Range, WindowData> entry : targetFile.getRanges().entrySet()) {
					if (mzRange.contains(entry.getKey().getMiddle())) {
						dutyCycleMap.put(entry.getKey(), entry.getValue());
					}
				}
				maskedFile.setRanges(dutyCycleMap);

			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		maskedFile.saveAsFile(outputFile);

		long stopTime = System.nanoTime();
		long duration = stopTime - startTime;
		System.out.println("maskByTargetWindows(): Time taken (ms) : " + duration / 1_000_000);

		return maskedFile;

	}

	// public static maskByTargetMz() {

	// }

}
