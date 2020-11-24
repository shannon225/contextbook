package edu.washington.gs.maccoss.encyclopedia;

import edu.washington.gs.maccoss.encyclopedia.utils.io.Version;

public enum ProgramType {
	EncyclopeDIA("EncyclopeDIA"), PecanPie("Walnut"), XCorDIA("XCorDIA"), CASiL("Thesaurus"), Global("Full EncyclopeDIA");
	
	private final String name;
	private ProgramType(String name) {
		this.name=name;
	}
	
	@Override
	public String toString() {
		return name;
	}
	
	public Version getVersion() {
		Version localVersion = getGlobalVersion();
		return localVersion;
	}

	public static Version getGlobalVersion() {
		String version = ProgramType.class.getPackage().getImplementationVersion();
		Version localVersion=new Version(version);
		return localVersion;
	}

	public String getCitation() {
		if (this==ProgramType.CASiL) {
			return "This is a <a href=http://villenlab.gs.washington.edu/>Villen Lab</a> and <a href=https://sites.google.com/a/uw.edu/maccoss/>MacCoss Lab</a> project from the University of Washington, <a href=http://www.gs.washington.edu/>Department of Genome Sciences</a>. "
					+ "You can cite this project as Searle BC et al. \"Thesaurus: quantifying phosphoprotein positional isomers.\" <a href=https://doi.org/10.1101/421214>Biorxiv 421214; doi: https://doi.org/10.1101/421214.</a> "
					+ "For more information please contact Brian Searle (bsearle@systemsbiology.org).";
		} else if (this==PecanPie) {
			return "This is a <a href=https://sites.google.com/a/uw.edu/maccoss/>MacCoss Lab</a> project from the University of Washington, <a href=http://www.gs.washington.edu/>Department of Genome Sciences</a>. "
					+ "This work is based on PECAN, originally published as Ting YS et al. \"PECAN: library-free peptide detection for data-independent acquisition tandem mass spectrometry data.\" <a href=https://doi.org/10.1038/nmeth.4390>Nat Methods. 2017;14:903-908.</a> "
					+ "You can cite this project as Searle BC et al. \"Chromatogram libraries improve peptide detection and quantification by data independent acquisition mass spectrometry.\" <a href=https://doi.org/10.1038/s41467-018-07454-w>Nat Commun. 2018;9:5128.</a> "
					+ "For more information please contact Brian Searle (bsearle@systemsbiology.org).";
		} else if (this==XCorDIA) {
			return "This is a <a href=https://sites.google.com/a/uw.edu/maccoss/>MacCoss Lab</a> project from the University of Washington, <a href=http://www.gs.washington.edu/>Department of Genome Sciences</a>. "
					+ "For more information please contact Brian Searle (bsearle@systemsbiology.org).";
		}
		return "This is a <a href=https://sites.google.com/a/uw.edu/maccoss/>MacCoss Lab</a> project from the University of Washington, <a href=http://www.gs.washington.edu/>Department of Genome Sciences</a>. "
			+ "You can cite this project as Searle BC et al. \"Chromatogram libraries improve peptide detection and quantification by data independent acquisition mass spectrometry.\" <a href=https://doi.org/10.1038/s41467-018-07454-w>Nat Commun. 2018;9:5128.</a> "
			+ "For more information please contact Brian Searle (bsearle@systemsbiology.org).";
	}

	public String getAboutMessage() {
		if (this==ProgramType.CASiL) {
			return "The human phosphorylation regulatory network represents a complex signaling cascade where proteins can be phosphorylated at multiple sites resulting in different functional states. Here we present Thesaurus, a hybrid search engine that detects and localizes novel positional isomers using site-specific fragment ions directly from data independent acquisition mass spectrometry experiments.";
		} else if (this==ProgramType.PecanPie) {
			return "<b>Walnut is just like PECAN, just with more wrinkles and slightly more bitter.";
		} else if (this==ProgramType.XCorDIA) {
			return "XCorDIA is XCorr for DIA. Duh.";
		} else {
			return "<b>en·cy·clo·pe·di·a</B>: a book or set of books giving information on many subjects or on many aspects of one subject and typically arranged alphabetically.";
		}
	}
}
