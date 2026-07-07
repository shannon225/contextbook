package edu.washington.gs.maccoss.encyclopedia.context;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class IsolationWindowReader {

	public static void main(String[] args) {
		if (args.length != 1) {
			System.err.println("Usage: java edu.washington.gs.maccoss.encyclopedia.context.IsolationWindowReader <massListFile>");
			System.exit(1);
		}
		ArrayList<IsolationWindow> windows = parseMassList(args[0]);
		System.out.println("Total windows read: " + windows.size());
	}

	/**
	 * Parse an assay / mass-list file into isolation windows.
	 *
	 * Two formats are accepted, auto-detected from the header:
	 *   1. Skyline / EncyclopeDIA isolation-list CSV (7 columns):
	 *        Compound, Formula, Adduct, m/z, z, RT Time (min), Window (min)
	 *   2. Targeted-assay TSV emitted by TargetedBootstrapper (8 columns):
	 *        compound, ..., m/z, charge, rtCenter, rtWindow, isDecoy
	 *
	 * The isDecoy flag defaults to false when the column is absent (format 1).
	 */
	public static ArrayList<IsolationWindow> parseMassList(String massListFile) {

		ArrayList<IsolationWindow> isolationWindows = new ArrayList<>();
		File massList = new File(massListFile);

		try (BufferedReader br = new BufferedReader(new FileReader(massList))) {

			String header = br.readLine();
			if (header == null) {
				throw new IOException("Mass list is empty: " + massListFile);
			}

			// Auto-detect delimiter: whichever character occurs more often in the header.
			int tabs   = countChar(header, '\t');
			int commas = countChar(header, ',');
			String delim = tabs >= commas ? "\t" : ",";

			String line;
			while ((line = br.readLine()) != null) {
				if (line.trim().isEmpty()) continue;

				String[] columns = line.split(delim, -1);
				if (columns.length < 7) {
					System.err.println("Skipping malformed mass-list row (expected >=7 fields, got "
							+ columns.length + "): " + line);
					continue;
				}

				String compound  = columns[0];
				double targetMz  = Double.parseDouble(columns[3]);
				byte   charge    = Byte.parseByte(columns[4]);
				float  rtCenter  = Float.parseFloat(columns[5]);
				float  rtWindow  = Float.parseFloat(columns[6]);
				boolean isDecoy  = (columns.length > 7) && Boolean.parseBoolean(columns[7]);

				float rtMin = (rtCenter - (rtWindow / 2)) * 60;
				float rtMax = (rtCenter + (rtWindow / 2)) * 60;

				isolationWindows.add(new IsolationWindow(compound, targetMz, charge, rtMin, rtMax, isDecoy));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		return isolationWindows;
	}

	private static int countChar(String s, char c) {
		int n = 0;
		for (int i = 0; i < s.length(); i++) if (s.charAt(i) == c) n++;
		return n;
	}
}
