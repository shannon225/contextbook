package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;

public abstract class QuantitativeSearchJobData extends SearchJobData {

	public QuantitativeSearchJobData(File diaFile, File featureFile, File outputFile, SearchParameters parameters, String version) {
		super(diaFile, featureFile, outputFile, parameters, version);
	}

	public abstract File getResultLibrary();
}