package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;

public abstract class DDASearchJobData extends AbstractSearchJobData {
	protected final StripeFileInterface diaFileReader;

	/**
	 * @param diaFileReader May be {@code null} if default generation behavior is desired
	 */
	public DDASearchJobData(File diaFile, StripeFileInterface diaFileReader, PercolatorExecutionData percolatorFiles, SearchParameters parameters, String version) {
		super(diaFile, percolatorFiles, parameters, version);

		this.diaFileReader=diaFileReader;
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