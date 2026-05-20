package edu.washington.gs.maccoss.encyclopedia.context;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class IsolationWindowReader {

	private static final String DELIM = "\t"; // How do I make this a csv?

	public static void main(String[] args) {

		String massListFile = "C:/Users/m334793/Documents/Library/targeted_bootstrapper_test/IT_100ngCurve_100p_masked0_assay.txt";

		ArrayList<IsolationWindow> windows = parseMassList(massListFile);
		
		System.out.println("Total windows read: " + windows.size());

	}

	// formatted as a mass list that is output when generating targeted assays with
	// encyclopedia
	public static ArrayList<IsolationWindow> parseMassList(String massListFile) {

		// Variables to fill in with the assay.csv entries
		ArrayList<IsolationWindow> isolationWindows = new ArrayList<>();
		File massList = new File(massListFile);

		try (BufferedReader br = new BufferedReader(new FileReader(massList))) {
	
			@SuppressWarnings("unused")
			String header = br.readLine();
			
			String line;
			while ((line = br.readLine()) != null) {
				String columns[] = line.split(DELIM, -1);
<<<<<<< HEAD
<<<<<<< HEAD
				String compound = columns[0];
=======
>>>>>>> f44678a1 (Added a class that will process features with Encyclopedia without running Percolator, and export them as pin.tsv files.)
=======
				String compound = columns[0];
>>>>>>> 8b8d896c (Added sequences to the mass lists for the TargetedBoostrapper class.)
				double targetMz = Double.parseDouble(columns[3]);
				byte charge = Byte.parseByte(columns[4]);
				float rtCenter = Float.parseFloat(columns[5]);
				float rtWindow = Float.parseFloat(columns[6]);
				boolean isDecoy = Boolean.parseBoolean(columns[7]);

				float rtMin = (rtCenter - (rtWindow / 2))*60;
				float rtMax = (rtCenter + (rtWindow / 2))*60;

//				boolean isDecoy = false;

				// Assemble each window
				IsolationWindow window = new IsolationWindow(compound, targetMz, charge, rtMin, rtMax, isDecoy);
				isolationWindows.add(window);
<<<<<<< HEAD
<<<<<<< HEAD
			System.out.println("Adding an mz at to the target list " + targetMz + " and RT " + rtCenter + " min " + rtMin/60 + " max " + rtMax/60 
=======
				System.out.println("Adding an mz at to the target list " + targetMz + " and RT " + rtCenter + " min " + rtMin/60 + " max " + rtMax/60 
>>>>>>> f44678a1 (Added a class that will process features with Encyclopedia without running Percolator, and export them as pin.tsv files.)
=======
//				System.out.println("Adding an mz at to the target list " + targetMz + " and RT " + rtCenter + " min " + rtMin/60 + " max " + rtMax/60 
>>>>>>> 8b8d896c (Added sequences to the mass lists for the TargetedBoostrapper class.)
//						+ "\nRTCenter: " + rtCenter
//						+ "\ntargetMz: " + targetMz
//						+ "\nrtStart: " + rtMinHmmm now
//						+ "\nrtStop: " + rtMax
//						);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		return isolationWindows;

	}
}
