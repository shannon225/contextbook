package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;

import java.io.File;

public abstract class QuantitativeSearchJobData extends AbstractSearchJobData {
	protected final StripeFileInterface diaFileReader;

	public QuantitativeSearchJobData(
			File diaFile,
			File featureFile,
			File outputFile,
			File decoyFile,
			SearchParameters parameters,
			String version
	) {
		this(diaFile, null, featureFile, outputFile, decoyFile, parameters, version);
	}

	/**
	 * @param diaFileReader May be {@code null} if default generation behavior is desired
	 */
	public QuantitativeSearchJobData(
			File diaFile,
			StripeFileInterface diaFileReader,
			File featureFile,
			File outputFile,
			File decoyFile,
			SearchParameters parameters,
			String version
	) {
		super(diaFile, featureFile, outputFile, decoyFile, parameters, version);

		this.diaFileReader = diaFileReader;
	}

	@Override
	public StripeFileInterface getDiaFileReader() {
		if (null != diaFileReader) {
			return diaFileReader;
		} else {
			return super.getDiaFileReader();
		}
	}

	public abstract File getResultLibrary();
}