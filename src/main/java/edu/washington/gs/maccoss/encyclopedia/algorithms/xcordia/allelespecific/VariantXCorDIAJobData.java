package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;

public class VariantXCorDIAJobData extends XCorDIAJobData {
	// GUI
	public VariantXCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, Optional<LibraryInterface> library, File diaFile, File fastaFile, XCorDIAOneScoringFactory taskFactory) {
		super(targetList, library, diaFile, fastaFile, taskFactory);
	}

	// command line
	public VariantXCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, Optional<LibraryInterface> library, File diaFile, File fastaFile, File outputFile, XCorDIAOneScoringFactory taskFactory) {
		super(targetList, library, diaFile, fastaFile, outputFile, taskFactory);
	}
	
	// internal
	public VariantXCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, Optional<LibraryInterface> library, File diaFile, StripeFileInterface diaFileReader, File fastaFile, PercolatorExecutionData percolatorFiles,
			XCorDIAOneScoringFactory taskFactory) {
		super(targetList, library, diaFile, diaFileReader, fastaFile, percolatorFiles, taskFactory);
	}

	public VariantXCorDIAJobData updateTaskFactory(XCorDIAOneScoringFactory taskFactory) {
		return new VariantXCorDIAJobData(getTargetList(), getLibrary(), getDiaFile(), getDiaFileReader(), getFastaFile(), getPercolatorFiles(), taskFactory);
	}

	public File getLocalizationFile() {
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile());
		return new File(absolutePath+".localizations.txt");
	}
}
